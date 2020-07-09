import requests
import datetime
import sys
import json 
import re
import subprocess
from requests.exceptions import Timeout,TooManyRedirects,RequestException
import config
import pandas as pd
from io import StringIO
import os.path
import logging
import jinja2

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

def getServerName(build="38"):
    server="http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server="http://rest.ensembl.org"
    return server

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

def makeOverlapVarQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/"
    return getServerName(build)+ext+"%s:%s-%s?feature=variation"  % (chrom,str(start),str(end))

def makeRefSeqQueryURL(chrom,start,end,build="38"):
    ext="/sequence/region/human/"
    return getServerName(build)+ext+"%s:%s..%s:1?"  % (chrom,str(start),str(end))

def makeRSQueryURL(rsID,build="38"):
    ext="/variant_recoder/homo_sapiens/"
    return getServerName(build)+ext+rsID+"?"

def makeHomologyURL(ID,species="mouse",build="38"):
    ext="/homology/id/%s?&target_species=%s&aligned=0&sequence=none&type=orthologues" %(ID,species)
    return getServerName(build)+ext

def makeGeneQueryURL(ID,build="38"):
    ext="/lookup/id/"
    return getServerName(build)+ext+ID

def makeGeneXQueryURL(ID,build="38"):
    ext="/xrefs/id/"
    return getServerName(build)+ext+"%s?all_levels=1" % ID

def makeGeneXQueryURL2(ID,build="38"):
    ext="/xrefs/id/"
    return getServerName(build)+ext+"%s?external_db=MGI" % ID

def makeRSPhenotypeQueryURL(build="38"):
    ext="/variation/homo_sapiens/"
    return getServerName(build)+ext+"?phenotypes=1"

def makePhenoOverlapQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/"
    return getServerName(build)+ext+"%s:%d-%d?feature=variation;variant_set=ph_variants" %(chrom,start,end)

def makeRsPhenotypeQuery2URL(rs,build="38"):
    ext="/variation/human/"
    return getServerName(build)+ext+"%s?pops=1;phenotypes=1" %rs

def makeRSListQueryURL(build="38"):
    ext="/variant_recoder/homo_sapiens/"
    return getServerName(build)+ext

def makeVepQueryURL(chrom,start,end,allele,strand="1",build="38"):
    ext="/vep/homo_sapiens/region/%s:%s-%s:1/%s?" % (chrom,str(start),str(end),allele)
    return getServerName(build)+ext

def makeGeneOverlapQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/%s:%d-%d?feature=gene" %(chrom,start,end)
    return getServerName(build)+ext

#------------------------------------------------------------------------------------------------------------------------------------

# VEP-style string
def getLocationString(chrom,pos,ref,alt):
    if len(ref)>1 and len(alt)==1:
        start=pos+1
        end=pos+len(ref)-1
        allele1=""
        allele2="-"
    elif len(alt)>1 and len(ref)==1:
        start=pos+1
        end=pos
        allele1=""
        allele2=alt[1:]
    elif len(ref)==1 and len(alt)==1:
        start=pos
        end=pos
        allele1=ref
        allele2=alt
    else:
        LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt),file=sys.stderr)
        return None

    return chrom+":"+str(start)+"-"+str(end)+":"+allele1+"/"+allele2

# given rsID, return a list of chr:pos
def rs2position(ID,build="38"):
    L=[]
    z=restQuery(makeRSQueryURL(ID,build=build))
    if z:
        for x in z:
            spdis=x["spdi"]
            for spdi in spdis:
                h=parseSPDI(spdi)
                p=h["pos"]
                c=h["chr"]
                z=next((x for x in L if x["chr"]==c and x["pos"]==p),None)
                if not z:
                    L.append({"chr":c,"pos":p})
    else:
        return None

    return L

#------------------------------------------------------------------------------------------------------------------------------------

# if alleles==True, also return alleles
def parseSPDI(string,alleles=False,build="38"):
    L=string.rsplit(":")
    c=L[0]
    m=re.search("NC_0+(\d+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=int(L[1])
    ref=L[2]
    alt=L[3]
    lref=len(ref)
    lalt=len(alt)
    isSNP=False
    base=None
    seq=None
    isNum=False

    # ref can be the length of the deleted sequence
    # or empty in case of simple insertions
    m=re.match("^(\d+)$",ref)
    if m:
        isNum=True
        if ref=="1" and lalt==1: # one base deleted, one inserted, i.e. variant is a SNP
            isSNP=True
    elif lref==1 and lalt==1:# SNP
        isSNP=True

    if alleles:
        # deleted sequence
        if isNum:
            delseq=getRefSeq(c,pos+1,pos+int(ref),build)
        else:
            delseq=ref

        if isSNP:
            return {"chr":c,"pos":pos+1,"ref":delseq[0],"alt":alt}
        else:
            # base before the deleted sequence
            base=getRefSeq(c,pos,pos,build)
            # deletion
            if len(delseq)>len(alt):
                if delseq.endswith(alt):
                    delseq2=delseq[0:len(delseq)-len(alt)]
                    return {"chr":c,"pos":pos,"ref":base+delseq2,"alt":base}
                else:
                    LOGGER.warning("parseSPDI: deletion: (%s %d %s %s) INS is not suffix of DEL" %(c,pos,ref,alt))
                    return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
            # insertion
            elif len(delseq)<len(alt):
                if lref==0:
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt}
                elif alt.endswith(delseq):
                    alt2=alt[0:len(alt)-len(ref)]
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt2}
                else:
                    LOGGER.warning("parseSPDI: insertion: (%s %d %s %s) DEL is not prefix of INS" %(c,pos,ref,alt))
                    return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
            # indel
            else:
                return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
    else:
        if isSNP:
            return {"chr":c,"pos":pos+1}
        else:
            return {"chr":c,"pos":pos}

def restQuery(URL,data=None,qtype="get",timeout=None):
    func=None

    if qtype=="get":
        func=requests.get
    elif qtype=="post":
        func=requests.post
    else:
        LOGGER.error("Query type %s has to be either \"get\" or \"post\"" %(qtype))
        return None

    r=None
    try:
        if qtype=="get":
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},timeout=timeout)
        else:
            if not data:
                LOGGER.error("POST query requires data")
                return None                
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},data=data,timeout=timeout)

        if not r.ok:
            LOGGER.error("Error %s occured (input: %s)" %(str(r.status_code),URL))
            return None

        try:
            ret=r.json()
            return ret
        except ValueError:
            LOGGER.error("JSON decoding error")
            return None

    except Timeout as ex:
        LOGGER.error("Timeout exception occured")
        return None
    except TooManyRedirects as ex:
        LOGGER.error("TooManyRedirects exception occured")
        return None
    except RequestException as ex:
        LOGGER.error("RequestException occured")
        return None

# check if provided input is a valid variant ID
# valid ID: either rsID or chr_pos_ref_alt
def checkID(id):
    m=re.search("^rs\d+",id)
    if m:
        return True
    
    m=re.search("^\d+_\d+_([ATGC]+)_([ATGC]+)",id)
    if m:
        if len(m.group(1))==1 and len(m.group(2))==1:
            return True

        if m.group(1).startswith(m.group(2)) or m.group(2).startswith(m.group(1)):
            return True
        else:
            return False
    else:
        return False

def chunks(L,n):
    """Yield successive n-sized chunks from L"""
    for i in range(0, len(L), n):
        yield L[i:i + n]

# for a given genomic region, return dataframe containing variants with phenotype annotations
def getVariantsWithPhenotypes(chrom,pos,window=config.WINDOW,build="38"):
    '''
    Input: chromosome, position, window (default: config.WINDOW), build (default: "38") 

    Output: pandas dataframe with the following columns:
    'SNPID', 'consequence', 'distance', 'phenotype', 'rsID', 'source'
    '''
    
    start=pos-window
    end=pos+window

    if start<1 : 
        start=1

    if end-start>5000000:
        LOGGER.error("Maximal region size allowed: 5Mbp")
        return None

    variants=restQuery(makePhenoOverlapQueryURL(chrom,start,end,build=build),qtype="get")
    #print(json.dumps(variants,indent=4,sort_keys=True))

    if not variants:
        return None

    if len(variants)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return None

    rsIDs = []
    for var in variants:
        rsIDs.append(var["id"])

    if len(rsIDs)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return None
    else:
        LOGGER.info("%d variant(s) with phenotypes were found in the region %s:%d-%d" %(len(rsIDs),chrom,start,end))

    output=[]
    i=0
    df = pd.DataFrame(columns=["ID","Consequence","Location","Phenotype","Source","Distance"])
    for L in chunks(rsIDs,config.BATCHSIZE):
        r=restQuery(makeRSPhenotypeQueryURL(build=build),data=list2string(L),qtype="post")
        if r:
            for rsID in r:
                for phenotype in r[rsID]["phenotypes"]:
                    m=re.match(".*phenotype\s+not\s+specified.*",phenotype["trait"])
                    if m:
                        continue
                    df.loc[i]=[rsID,r[rsID]["most_severe_consequence"].replace("_"," "),r[rsID]["mappings"][0]["seq_region_name"]+":"+str(r[rsID]["mappings"][0]["start"]),phenotype["trait"],phenotype["source"],str(abs(pos - int(r[rsID]["mappings"][0]["start"])))]
                    i+=1
            # # TODO: check all possible sources
            # # Generate phenotype link based on source:
            # if phenotype["source"] == "ClinVar":
            #     link = "https://www.ncbi.nlm.nih.gov/clinvar/?term="+rsID
            # elif phenotype["source"] == "HGMD-PUBLIC":
            #     link = "http://www.hgmd.cf.ac.uk/ac/gene.php?gene="+phenotype["genes"]
            # elif "NHGRI-EBI" in phenotype["source"]:
            #     link = "https://www.ebi.ac.uk/gwas/search?query="+rsID
            # elif phenotype["source"] == "OMIM":
            #     link = "https://www.omim.org/entry/"+phenotype["study"][4:]
            # else:
            #     link = "http://"+buildstr+"ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;v=CM106680;vdb=variation"

            # TODO: check key availability

    return df

# ===========================================================================================================================

# get general information about a variant, given rsID:
# 
# MAF, minor allele, variant class, most severe consequence
# 1KG phase 3 population allele frequencies
# mappings: chr, pos, ref, alt
# phenotype data: trait, source, risk allele
# variant's clinical significance
# TODO: remove MAF and "minor_allele" ?
def getVariantInfo(rs,build="38"):
    res=dict()

#------------------- general information ---------------

    data=restQuery(makeRsPhenotypeQuery2URL(rs,build))
    #print(json.dumps(data,indent=4,sort_keys=True))

    if not data:
        return res

    res["minor_allele"]=data["minor_allele"]
    res["MAF"]=data["MAF"]
    res["rsID"]=rs
    res["class"]=data["var_class"]
    res["consequence"]=data["most_severe_consequence"]
    if "synonyms" in data:
        res["synonyms"]=data["synonyms"]
    else:
        res["synonyms"]=[]

#------------------- mappings----------------------

    mappings=list()

    z=restQuery(makeRSQueryURL(rs,build=build))
    for x in z:
        spdis=x["spdi"]
        for spdi in spdis:
            h=parseSPDI(spdi,alleles=True)
            ref=h["ref"]
            alt=h["alt"]
            p=h["pos"]
            c=h["chr"]
            mappings.append({"chr":c,"pos":p,"ref":ref,"alt":alt,"sift_score":"NA","sift_prediction":"NA","polyphen_score":"NA","polyphen_prediction":"NA"})

#------------------ population data ----------------

    population_data=list()

    for pop in data["populations"]:
        pop_name = pop["population"].split(":")
        if pop_name[0] == "1000GENOMES" and pop_name[1] == "phase_3":
            name=pop_name[2]
            try:
                z=next(x for x in population_data if name==x["population"])
                z["frequency"][pop["allele"]]=pop["frequency"]
            except:
                population_data.append({"population":name,"frequency":{pop["allele"]:pop["frequency"]}})

#------------------ phenotype data -------------------

    phenotype_data=list()

    for p in data["phenotypes"]:
        trait=p["trait"] if "trait" in p else "NA"
        source=p["source"] if "source" in p else "NA"
        risk=p["risk_allele"] if "risk_allele" in p else "NA"
        if trait:
            phenotype_data.append({"trait":trait,"source":source,"risk_allele":risk})


#------------------ clinical significance -------------------

    clinical_significance=list()

    if "clinical_significance" in data:
        for cs in data["clinical_significance"]:
            if cs != "other" and cs != "not provided":
                clinical_significance.append(cs)

#-----------------------------------------------------

    res["mappings"]=mappings
    res["population_data"]=population_data
    res["phenotype_data"]=phenotype_data
    res["clinical_significance"]=clinical_significance

    return res

# ===========================================================================================================================

# return a part of the reference genome
def getRefSeq(chrom,start,end,build="38"):
    r=restQuery(makeRefSeqQueryURL(chrom,start,end,build))
    if r:
        return r["seq"]
    else:
        return None

# ===========================================================================================================================

# given variant ID, try to get variant type
# also check if reference sequence and alleles are not in contradiction with each other
def getVariantType(varid,build="38"):
    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
    chrom=m.group(1)
    pos=int(m.group(2))
    a1=m.group(3)
    a2=m.group(4)

    vartype=None
    seq=getRefSeq(chrom,pos,pos+max([len(a1),len(a2)])-1,build=build)
    if len(a1)==1 and len(a2)==1:
        if seq==a1 or seq==a2:
            vartype="SNP"
        else:
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: SNP alleles don't match reference sequence", file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Variant: %s" % varid, file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Reference: %s\n" % seq, file=sys.stderr)
            sys.stderr.flush()
            return None
    elif len(a1)<len(a2):
        if not seq.startswith(a1):
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: variant alleles don't match reference sequence", file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Allele 1  : %s" % a1, file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Allele 2  : %s" % a2, file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Reference : %s\n" % seq, file=sys.stderr)
            sys.stderr.flush()
            return None
        if a2!=seq:
            vartype="INS"
    else:
        if not seq.startswith(a2):
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: variant alleles don't match reference sequence", file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Allele 1  : %s" % a1, file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Allele 2  : %s" % a2, file=sys.stderr)
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantType: Reference : %s\n" % seq, file=sys.stderr)
            sys.stderr.flush()
            return None
        if a1!=seq:
            vartype="INS"

    return vartype

# ===========================================================================================================================

def stringdiff(s1,s2):
    L=[]

    t1=s2
    t2=s1

    if len(s1)<len(s2):
        t1=s1
        t2=s2

    for i in range(len(t1)+1):
        x1=t1[i:len(t1)]
        x2=t1[0:i]
        if t2.startswith(x2) and t2.endswith(x1):
            L.append(t2[i:len(t2)-(len(t1)-i)])

    return L

# ===========================================================================================================================

# returns a set of matching rsIDs for a given variant ID
def id2rs(varid,build="38"):
    if varid.startswith("rs"):
        return varid

    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
    chrom=m.group(1)
    pos=int(m.group(2))
    a1=m.group(3)
    a2=m.group(4)

    # in case of indels, pull all variants around the given variant
    window=5

    batchsize=100

    S=set()
    if len(a1)==1 and len(a2)==1: # SNP
        r=restQuery(makeOverlapVarQueryURL(chrom,pos,pos,build=build))
        for v in r:
            if a1 in v["alleles"] and a2 in v["alleles"]:
                S.add(v["id"])
    else:
        # difference between a1 and a2
        # seq0=""
        # if len(a1)>len(a2):
        #     seq0=a1[len(a2):len(a1)]
        # elif len(a1)<len(a2):
        #     seq0=a2[len(a1):len(a2)]

        r=restQuery(makeOverlapVarQueryURL(chrom,pos-window,pos+window,build=build))
        if not r:
            return S

        for v in r:
            z=restQuery(makeRSQueryURL(v["id"],build=build))
            if not z:
                continue

            for x in z:
                spdis=x["spdi"]
                var=x["id"][0]
                for spdi in spdis:
                    h=parseSPDI(spdi,alleles=True)
                    ref=h["ref"]
                    alt=h["alt"]
                    p=h["pos"]
                    c=h["chr"]

                    if p!=pos:
                        continue

                    if len(ref)==1 and len(alt)==1:
                        continue

                    if (ref==a1 and alt==a2) or (ref==a2 and alt==a1):
                    #if seq0 in stringdiff(alt[1:len(alt)],ref[1:len(ref)]):
                        S.add(var)
                        break
    return S

# ===================================================== GTEx RELATED STUFF ============================================

# given gene ID (variant ID), retreive all variant (gene) data associated with the gene (variant): tissue, p-value
def parseGTEx(chrom,start,end,ID):
    '''
    Input  : gene ID or variant ID (format: chr_pos_ref_alt)
    Output : dictionary with variant IDs (gene IDs) as keys and dictionaries with "tissue", "p-value", "beta", "SE", "dist" as values
    '''
    filename=config.GTEX_BED
    query = "tabix %s %s:%d-%d | awk -v v=%s \'BEGIN{FS=\"\t\";}$4==v{print $0;}\'" %(filename,chrom,start,end,ID)
    output = subprocess.Popen(query.strip(),universal_newlines=True,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    D=dict()
    for line in output.stdout.readlines():
        fields=line.strip().split("\t")
        data=fields[4].split(":")
        ID2=data[0]
        tissue=data[1]
        pval=str("%.4e" % float(data[2]))
        beta=str(round(float(data[3]),4))
        SE=str(round(float(data[4]),4))
        dist=data[5]
        if not ID2 in D:
            D[ID2]=[]
        D[ID2].append({"tissue" : tissue,"p-value" : pval,"beta" : beta,"SE" : SE,"dist" : dist})

    if len(D)==0:
        LOGGER.info("No eQTL signals were found for %s" %(ID))
    
    return D

# ===================================================== UNIPROT RELATED STUFF ============================================

def getUniprotData(ID):
    URL = config.UNIPROT_URL
    URL += "?query=id:" + ID
    URL += "&columns=id%2Ccomment%28FUNCTION%29%2C" # Function
    URL += "comment%28SUBUNIT%29%2C" # Subunit
    URL += "comment%28DEVELOPMENTAL%20STAGE%29%2C" # Developmental stage
    URL += "comment%28TISSUE%20SPECIFICITY%29%2C" # Tissue specificity
    URL += "comment%28CATALYTIC%20ACTIVITY%29%2C" # Catalytic activity
    URL += "comment%28DISRUPTION%20PHENOTYPE%29%2C" # Disruption phenotype
    URL += "comment%28SUBCELLULAR%20LOCATION%29%2C" # Subcellular localization
    URL += "comment%28DISEASE%29%2C" # Disease
    URL += "entry%20name&format=tab" # tab delimited format returned

    r=requests.get(URL)

    data=pd.read_csv(StringIO(r.content.decode("utf-8")),sep="\t",header=0,keep_default_na=False)

    UniprotData = {
        "Name" : data["Entry name"][0],
        "Disease" : data["Involvement in disease"][0],
        "Function": data["Function [CC]"][0],
        "Entry" : data["Entry"][0],
        "Subunit" : data["Subunit structure [CC]"][0],
        "Phenotype" : data["Disruption phenotype"][0],
        "Location" : data["Subcellular location [CC]"][0],
        "Tissue" : data["Tissue specificity"][0],
        "Development" : data["Developmental stage"][0]
    }

    return UniprotData

# ================================================ GWAS CATALOG ===========================================================

def gene2gwas(gene_name):
    '''
    This function performs a query for gwas signals in a local gwas datafile.
    All reported and mapped genes will give a hit.

    Output: list of dictionaries
    '''

    gwas_file = config.GWAS_FILE
    query = "zcat %s | fgrep -iw %s" % (gwas_file,gene_name)

    gwas_data = []
    output = subprocess.Popen(query,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    for line in output.stdout.readlines():
        fields = line.split("\t")

        gwas_data.append({
            "rsID": fields[21],
            "p-value": fields[27],
            "Consequence": fields[24],
            "Reported genes": fields[13],
            "Mapped genes": fields[14],
            "trait": fields[7],
            "PMID": fields[1],
            "URL": fields[5],
            "Author": fields[2]
        })

    if len(gwas_data)==0:
        LOGGER.info("No GWAS signas associated with %s have been found" %(gene_name))

    return gwas_data


def getGwasHits(chrom,pos,window=500000):
    '''
    Retrives a list of all gwas signals within a window around a given position

    Input: chromosome, position, window size (default: 500000)

    Output: [ {"rsID","SNPID","trait","p-value","PMID","distance"}]
    '''

    gwas_file=config.GWAS_FILE_VAR

    if not os.path.isfile(gwas_file):
        LOGGER.error("GWAS catalog file (%s) not found" % gwas_file)
        return None

    start=pos-window
    if start<1:
        start=1
    end=pos+window

    L=[]
    df=pd.read_table(gwas_file,sep="\t",header=0,compression="gzip")
    for index, row in df[(df["CHR_ID"]==chrom) & (df["CHR_POS"].astype(int)>start) & (df["CHR_POS"].astype(int)<end)].iterrows():
        rsID=row["SNPS"]
        snpid=row["CHR_ID"]+"_"+str(row["CHR_POS"])
        trait=row["DISEASE/TRAIT"]
        pval=row["P-VALUE"]
        pmid=row["PUBMEDID"]
        dist=int(row["CHR_POS"])-pos
        if abs(dist)>1000:
            dist=str(dist//1000)+"kbp"
        else:
            dist=str(dist)+"bp"
        L.append({"rsID":rsID,"SNPID":snpid,"trait":trait,"p-value":pval,"PMID":pmid,"distance":dist})

    return L

# ==============================================================================================================================

# def get_GWAVA_score(variant_data):
#     '''
#     This function calculates the gerp and gwava scores for a given variation.
#     Of course it calls the gwava on farm, and parses the result.
#     '''

#     # Using variant_data for input:
#     bed_string = "chr%s\t%s\t%s\t%s\n" % (variant_data["chromosome"], int(variant_data["start"]) - 1, variant_data["start"],variant_data["rsID"])

#     # temporary input filename:
#     filename = "/tmp/chr%s.bed" % variant_data["location"].replace(":", "_")

#     # temporary output filename with gwava annotation:
#     annot_filename = filename+"_ann"

#     # temporary output filename with gwava prediction:
#     gwava_filename = filename+"_gwava"

#     # Saving primary input file:
#     f = open( filename, 'w')
#     f.write(bed_string)
#     f.close()

#     # now we have to run gwava:
#     GWAVA_dir = config.GWAVA_DIR

#     #query = "bash -O extglob -c \'python %s/src/gwava_annotate.py %s %s\'" %(GWAVA_dir, filename, annot_filename)
#     query = "python %s/src/gwava_annotate.py %s %s" %(GWAVA_dir, filename, annot_filename)

#     PATH = "/software/hgi/pkglocal/samtools-1.2/bin:/software/hgi/pkglocal/vcftools-0.1.11/bin:/software/hgi/pkglocal/tabix-git-1ae158a/bin:/software/hgi/pkglocal/bcftools-1.2/bin:/nfs/team144/software/ensembl-releases/75/ensembl-tools/scripts/variant_effect_predictor:/nfs/team144/software/bedtools2/bin:/nfs/team144/software/scripts:/nfs/users/nfs_d/ds26/bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/software/bin"

#     # Submit query:
#     output = subprocess.Popen(query.strip(),
#                               shell=True,
#                               stdout=subprocess.PIPE,
#                               stderr=subprocess.STDOUT,
#                               env={'GWAVA_DIR': GWAVA_dir,
#                                    "PYTHONPATH": "/nfs/team144/software/anaconda/lib/python2.7/site-packages",
#                                    "PATH": PATH}).wait()

#     # Once the annotation run is completed, let's run the prediction:
#     query = "python %s/src/gwava.py tss %s %s" %(GWAVA_dir, annot_filename, gwava_filename)

#     # Submit query:
#     output = subprocess.Popen(query.strip(),
#                               shell=True,
#                               stdout=subprocess.PIPE,
#                               stderr=subprocess.STDOUT,
#                               env={'GWAVA_DIR': GWAVA_dir,
#                                    'PYTHONPATH': '/nfs/team144/software/anaconda/lib/python2.7/site-packages',
#                                    'PATH': PATH}).wait()

#     # Once the queries are returned, we have to parse the output:
#     # From the annotation file, we retrive 147 - avg_gerp, 148 - gerp

#     # Reading annotation file:
#     for line in open(annot_filename, 'r'):
#         if "start" in line: continue

#         avg_gerp = line.split(",")[147]
#         gerp = line.split(",")[148]

#     # Reading prediction file:
#     for line in open(gwava_filename, 'r'):
#         gwava = line.strip().split("\t")[4]

#     # Updating variant
#     variant_data["avg_gerp"] = round(float(avg_gerp), 4)
#     variant_data["gerp"] = gerp
#     variant_data["gwava"] = gwava

#     # removing temporary files:
#     os.remove(annot_filename)
#     os.remove(filename)
#     os.remove(gwava_filename)

#     return variant_data

# ==============================================================================================================================

# Download gene information from the Ensembl:
def getGeneInfo (ID,build="38"):
    '''
    This function retrieves gene related information

    Input: Ensembl stable ID
    Output: dictionary with retrieved information
    '''
    response=restQuery(makeGeneQueryURL(ID,build=build))

    gene_info=dict()
    gene_info["source"] = response["source"]
    gene_info["id"] = response["id"]
    gene_info["start"] = response["start"]
    gene_info["end"] = response["end"]
    gene_info["assembly_name"] = response["assembly_name"]
    try:
        gene_info["description"] = response["description"].split("[")[0]
    except:
        gene_info["description"] = "NA"
    gene_info["name"] = response["display_name"]
    gene_info["type"] = response["biotype"]
    gene_info["strand"] = response["strand"]
    gene_info["chromosome"] = response["seq_region_name"]

    return gene_info

# ==============================================================================================================================

# Download cross references to a gene based on Ensembl annotation:
def getGeneXrefs (ID,build="38"):
    '''
    This function retrieves cross-references from Ensembl
    Input  : gene ID, build(default: "38")
    Output : dictionary with keys "MIM disease","MIM gene","GO","GOSlim GOA","UniProtKB/Swiss-Prot","Human Protein Atlas","ChEMBL"
                        and lists of tuples ("primary_id","description") as values  
    '''

    r=restQuery(makeGeneXQueryURL(ID,build=build))
    if not r:
        LOGGER.info("No cross references found for %s" %(ID))
        return  None

    xrefs = {
        "MIM disease" : [],
        "MIM gene" : [],
        "GO" : [],
        "GOSlim GOA" : [],
        "UniProtKB/Swiss-Prot" : [],
        "Human Protein Atlas" : [],
        "ChEMBL" : [],
    }

    for xref in r:        
        db=xref['db_display_name']
        if db in xrefs:
            x=(xref["primary_id"],xref["description"])
            if x not in xrefs[db]:
                xrefs[db].append(x)

    return xrefs

# ==============================================================================================================================

# Retrieve data from the OMIM database
# def getOmimInfo (cross_refs):
#     '''
#     Required information:
#         - list of OMIM ID-s

#     Retrieved information:
#         - OMIM entry name
#         - text field name
#         - text field (that includes references as well)

#     Retrieved data structure (one dict for all queried entries):
#     dict{
#         OMIM_ID :{
#             title : OMIM_preferredTitle
#             text : {
#                 OMIM_textSectionTitle : OMIM_textSectionContent
#     }}}

#     More information how the OMIM API works see: http://omim.org/help/api
#     The API key will expire in 2016.11.11, at that time a new key have to be required.
#     '''

#     # Extracting OMIM ID from the cross refereces data:
#     MIM_IDs = [x[1] for x in cross_refs["MIM disease"]]
#     MIM_IDs.extend([x[1] for x in cross_refs["MIM gene"]])

#     # The function returns at this point if there is not MIM IDs in the crossref
#     if len(MIM_IDs) == 0:
#         return "No OMIM entry for this gene."

#     # Constructing query string:
#     URL = config.OMIM_URL

#     for ID in MIM_IDs:
#         URL += 'mimNumber=%s&' % ID

#     # Retrieving the following fields:
#     URL += 'include=text&'
#     URL += 'include=allelicVariantList&'
#     URL += 'include=referenceList&'

#     # Adding API key:
#     URL += config.OMIM_APIKEY

#     # Required format is pyton (although it is a string indeed)ormat=python
#     URL += "format=python&"

#     # Downloading page:
#     page = requests.get(URL)

#     # Reconstructingh dictionary from the returned string:
#     OMIM_entry  = ast.literal_eval(page.content)

#     # hash to fill in:
#     OMIM_data = {}

#     # Parsing hash:
#     for entry in OMIM_entry['omim']['entryList']:
#         # GEt OMIM ID:
#         ID = entry["entry"]["mimNumber"]
#         OMIM_data[ID] = {}

#         # Get OMIM name:
#         OMIM_data[ID]['title'] = entry["entry"]["titles"]['preferredTitle']

#         # Get OMIM text:
#         OMIM_data[ID]['text'] = {}
#         for fields in entry['entry']["textSectionList"]:
#             OMIM_data[ID]['text'][fields['textSection']['textSectionTitle']] = fields['textSection']['textSectionContent']

#         # now we have to parse allelic variants:
#         # print stuff['omim']['entryList'][0]['entry']['allelicVariantList'][0]['allelicVariant'].keys()
#         # ['status', 'name', 'dbSnps', 'text', 'mutations', 'number', 'alternativeNames', 'clinvarAccessions']

#         try:
#             OMIM_data[ID]['variations'] = {}
#             for variations in entry['entry']['allelicVariantList']:
#                 OMIM_data[ID]['variations'][variations['allelicVariant']['dbSnps']] = variations['allelicVariant']['text']
#         except:
#             continue

#     return OMIM_data

# ==============================================================================================================================

# Retrievig data from the EBI's gene expression atlas:
# def getGxaData(ID):
#     '''
#     This function downloads information from the EBI's gene expression atlas.
#     Input: Ensembl ID
#     Output: pandas core series with the 10 most highly expressing tissues.
#     '''
#     # the URL pointing to the gene expression atlas:
#     URL = config.GXA_URL % ID

#     # The downloaded data will be read as a pandas dataframe:
#     try:
#         df = pd.read_csv(URL, sep='\t', comment="#")
#     except:
#         return ("No gene expression information!", "No gene information.")

#     cleanExperiments = []
#     for exp in df.Experiment.tolist():
#         m = re.search('Tissues -\s+\d*(.+)', exp)
#         cleanExperiments.append(m.groups()[0].strip())

#     df.Experiment = cleanExperiments

#     IndexNames = df.Experiment.tolist()

#     # These are the main sources that we focus on:
#     sources = ["GTEx", "FANTOM5 project - adult", "FANTOM5 project - fetal"]
#     Levels = {}

#     # Extracting top 10 tissues for each sources:
#     for source in sources:
#         try:
#             fieldName = filter(lambda x : source in x, IndexNames)[0]
#             Index = IndexNames.index(fieldName)

#             # Excising the GTEx data from dataframe:
#             GTEx = df.ix[Index]

#             GTEx = GTEx.sort_values(ascending=False, na_position='last')[0:11]
#             Levels[fieldName] = []
#             for tissue, value in GTEx.iteritems():
#                 try:
#                     if float(value) > 0 : Levels[fieldName].append([tissue, value])
#                 except:
#                     pass
#         except:
#             pass

#     return (Levels, df)

# ========================================================== MOUSE STUFF ===========================================================

def getMouseID(human_ID,build="38"):
    '''Looking up mouse gene ID of a given human gene ID'''

    data=restQuery(makeHomologyURL(human_ID,build=build,species="mouse"))
    #print(json.dumps(data,indent=4,sort_keys=True))

    mouse_IDs = {}
    if len(data["data"])==0 or len(data["data"][0]["homologies"])==0:
        LOGGER.info("No mouse cross-reference for %s" %(human_ID))
        return mouse_IDs

    for homolog in data["data"][0]["homologies"]:
        try:
            z=restQuery(makeGeneQueryURL(homolog["target"]["id"],build=build))
            #print(json.dumps(z,indent=4,sort_keys=True))

            d=""
            m=re.search("(.*)\s+\[.*\]",z["description"])
            if m:
                d=m.group(1)
            mouse_IDs[z["id"]] = [d,z["display_name"]]
        except:
            continue

    return mouse_IDs

def getMgiID(mouse_ID,build="38"):
    '''Looking up MGI cross reference for a given mouse gene ID'''

    data=restQuery(makeGeneXQueryURL2(mouse_ID,build=build))
    #print(json.dumps(data,indent=4,sort_keys=True))

    # Now from the returned data we have to pick all possible homolgue IDs:
    for d in data:
        try:
            return d["primary_id"]
        except:
            continue

    LOGGER.info("MGI ID for %s not found" % (mouse_ID))
    return ""

def getMgiPhenotypes(MGI_ID):
    '''Returns phenotype information stored on http://www.informatics.jax.org/ '''

    URL = "http://www.informatics.jax.org/allele/report.txt?markerId=%s" % MGI_ID

    # Return all associated allele information:
    try:
        df=pd.read_csv(URL, sep="\t")

        # Drop unnecessary columns:
        df.drop([ "Allele Symbol", "Chromosome",  "Synonyms","Allele Attributes", "Transmission", "Unnamed: 10"], axis=1, inplace=True)
        df.columns = ["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease"]
        return df
    except:
        LOGGER.info("No phenotype was found for %s" %(MGI_ID))
        return pd.DataFrame(columns=["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease"])

def getMousePhenotypes(geneID,build="38"):
    '''returning mouse phenotype given human gene ID'''
    # Returning all mouse homologue IDs:
    mouse_gene_IDs=getMouseID(geneID,build=build)

    full_dataframe=pd.DataFrame(columns=["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease","mouse gene ID","MGI ID","mouse gene name","mouse gene description"])

    if len(mouse_gene_IDs)==0:
        return full_dataframe

    MGI_IDs=dict()
    for mouse_gene_ID in mouse_gene_IDs:
        MGI_IDs[mouse_gene_ID] = getMgiID(mouse_gene_ID,build=build)

    # Once we have all the MGI identifiers, we retrieve all the phenotypes:
    for mouse_id, mgi_id in MGI_IDs.items():
        df=getMgiPhenotypes(mgi_id)

        # Adding extra columns for the record:
        df["mouse gene ID"] = mouse_id
        df["MGI ID"] = mgi_id
        df["mouse gene name"] = mouse_gene_IDs[mouse_id][1]
        df["mouse gene description"] = mouse_gene_IDs[mouse_id][0]
        full_dataframe = pd.merge(full_dataframe, df, how="outer")

    return full_dataframe

# ======================================================= APPRIS ===========================================================

def getApprisInfo(gene_ID):
    '''
    This function returns transcript information from the APPRIS webserver based on Ensembl gene ID
    The output is a dictionary with detailed information of all annotated transcripts

    Link: http://appris.bioinfo.cnio.es/
    '''

    #LOGGER.info("Calling getApprisInfo")

    # hardcoded assembly
    URL = "http://apprisws.bioinfo.cnio.es:80/rest/exporter/id/homo_sapiens/%s?format=json&db=b38" % gene_ID;

    data=dict()
    r = requests.get(URL)
    if not r.ok:
        LOGGER.info("Query failed for gene: %s (URL: %s)" % (gene_ID,URL))
        return data

    decoded = r.json()
    #print(json.dumps(decoded,indent=4,sort_keys=True))
    for transcript in decoded:
        tid=transcript["transcript_id"]
        if not tid in data:
            data[tid] = {"reliability" : "NA","name" : "NA","length_aa" : "NA","length_na" : "NA","type" : []}
        if "reliability" in transcript:
            data[tid]["reliability"]=transcript["reliability"]
        if "name" in transcript:
            data[tid]["name"]=transcript["name"]
        if "length_aa" in transcript:
            data[tid]["length_aa"]=transcript["length_aa"]
        if "length_na" in transcript:
            data[tid]["length_na"]=transcript["length_na"]
        if "type" in transcript:
            t=transcript["type"]
            if t not in data[tid]["type"]:
                data[tid]["type"].append(t)

    return data

# ======================================================= VEP ===========================================================

def getVepDF(mappings):
    '''
    For a given list of variant mappings (containing chr/pos/ref/alt information), 
    return two merged dataframes
    
    Input  : list of mappings
    Output : dictionary ("regulatory": regulatory DataFrame,"transcript": transcript DataFrame)
    '''

    LOGGER.debug("Input: %d mappings" %  len(mappings))

    output=dict()
    tr_data=pd.DataFrame(columns=["Gene ID","Transcript ID","Impact","Consequence","Isoform"])
    reg_data=pd.DataFrame(columns=["Biotype","Regulatory feature ID","Impact","Consequence"])
    for m in mappings:
        data=getVepData(m)
        tr_data=pd.concat([tr_data,vepTranscript2df(data)]).drop_duplicates().reset_index(drop=True)
        reg_data=pd.concat([reg_data,vepRegulatory2df(data)]).drop_duplicates().reset_index(drop=True)
        
    output["transcript"]=tr_data
    output["regulatory"]=reg_data
    
    return output

# Retrives data from the variant effect predictor server
# Modifies uinput mapping_data by adding SIFT/PolyPhen scores and predictions, if available
def getVepData(mapping_data):
    '''
    This function returns variant's predicted effect based on chromosome, position and the alternate allele.

    Input: variant data, a dictionary with chr, pos, ref, alt

    Output: dictionary with keys "transcript", "regulatory"
    "transcript" : list of dictionaries with keys "ID", "impact", "gene_symbol", "gene_id", "consequence", "principal"
    "regulatory" : list of dictionaries with keys "impact", "biotype", "ID", "consequence"    
    '''

    chrom = mapping_data["chr"]
    pos = mapping_data["pos"]
    ref = mapping_data["ref"]
    alt = mapping_data["alt"]

    start=-1
    end=-1
    allele=""

    if len(ref)>1 and len(alt)==1:
        start=pos+1
        end=pos+len(ref)-1
        allele="-"
    elif len(alt)>1 and len(ref)==1:
        start=pos+1
        end=pos
        allele=alt[1:]
    elif len(ref)==1 and len(alt)==1:
        start=pos
        end=pos
        allele=alt
    else:
        LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt),file=sys.stderr)
        return None

    VEP=restQuery(makeVepQueryURL(chrom,start,end,allele))
    #return data

    sift_score=None
    polyphen_score=None
    sift_prediction=""
    polyphen_prediction=""

    VEP_data=dict()
    VEP_data["transcript"]=[]
    VEP_data["regulatory"]=[]

# ------------------------------------------------------------------------------------------------------------------

    if "transcript_consequences" in VEP[0]:
        for t in VEP[0]["transcript_consequences"]:
            VEP_data["transcript"].append({"ID":t["transcript_id"],"impact":t["impact"],"gene_symbol":t["gene_symbol"],"gene_id":t["gene_id"],"consequence":t["consequence_terms"][0],"principal":"NA"})

            if "polyphen_score" in t:
                if polyphen_score:
                    if t["polyphen_score"]>polyphen_score:
                        polyphen_score=t["polyphen_score"]
                        polyphen_prediction=t["polyphen_prediction"]
                else:
                    polyphen_score=t["polyphen_score"]
                    polyphen_prediction=t["polyphen_prediction"]

            if "sift_score" in t:
                if sift_score:
                    if t["sift_score"]<sift_score:
                        sift_score=t["sift_score"]
                        sift_prediction=t["sift_prediction"]
                else:
                    sift_score=t["sift_score"]
                    sift_prediction=t["sift_prediction"]

# ------------------------------------------------------------------------------------------------------------------

    if sift_score is not None:
        mapping_data["sift_score"]=sift_score
        mapping_data["sift_prediction"]=sift_prediction
        
    if polyphen_score is not None:
        mapping_data["polyphen_score"]=polyphen_score
        mapping_data["polyphen_prediction"]=polyphen_prediction
        
# ------------------------------------------------------------------------------------------------------------------

    if "regulatory_feature_consequences" in VEP[0]:
        for r in VEP[0]['regulatory_feature_consequences']:
            VEP_data["regulatory"].append({"impact":r["impact"],"biotype":r["biotype"],"ID":r["regulatory_feature_id"],"consequence":r["consequence_terms"]})

# ------------------------------------------------------------------------------------------------------------------

    all_genes=[]
    for t in VEP_data["transcript"]:
        if not t["gene_id"] in all_genes:
            all_genes.append(t["gene_id"])

    for gene_id in all_genes:
        appris_data=getApprisInfo(gene_id)

        for t in VEP_data["transcript"]:
            if t["gene_id"]==gene_id:
                if t["ID"] in appris_data:
                    if "reliability" in appris_data[t["ID"]]:
                        t["principal"]=appris_data[t["ID"]]["reliability"]

    return VEP_data

# ======================================================= PUBMED ===========================================================

# searches for publications containing given rsID and its synonyms, returns merged dataframe
def getPubmedDF(rsID, synonyms):
    '''
    Searches for publications containing given rsID and its synonyms, returns merged dataframe
    Input: rsID, list of synonyms
    Output: dataframe with columns "firstAuthor", "journal", "year", "URL", "title"

    '''

    LOGGER.debug("Creating DF for %s" % rsID)
    df=pubmed2df(getPubmed(rsID))
    if len(synonyms)==0:
        return df
    else:
        for v in synonyms:
            LOGGER.debug("Creating DF for %s" % v)
            df2=pubmed2df(getPubmed(v))
            LOGGER.debug("Merging")
            df=pd.concat([df,df2]).drop_duplicates().reset_index(drop=True)

    return df

def getPubmed(rsID):
    '''
    This function returns a list of PMIDs of those publications where the given rsID was mentioned
    Up to 1000 IDs are returned

    Input  : rsID
    Output : dictionary with PMID as keys and dictionaries {"firstAuthor", "title", "journal", "year", "URL"} as values
    '''
    decoded=restQuery(config.PUBMED_URL_VAR % (rsID))
    #json.dumps(decoded,indent=4,sort_keys=True)
    pubmed_IDs=decoded["esearchresult"]["idlist"]

    publication_data = {}
    for ID in pubmed_IDs:
        r=requests.get(config.PUBMED_URL_PMID % (ID))
        decoded=r.json()

        # Extracting data:
        publication_data[ID] = {
            "firstAuthor" : decoded["result"][ID]['sortfirstauthor'],
            "title" : decoded["result"][ID]['title'],
            "journal" : decoded["result"][ID]['fulljournalname'],
            "year" : decoded["result"][ID]['epubdate'].split(" ")[0],
            "URL" : "http://www.ncbi.nlm.nih.gov/pubmed/"+ID,
        }

    return publication_data

# ======================================================================================================================

# Based on a genomic location, this function retrieves a list of genes within a window around it
def getGeneList(chrom,pos,window=1000000,build="38"):
    '''
    Based on the submitted chromosome and position, this function returns
    all genes within a window

    Input: chromosome, position, window (default: 1Mbp), build (default: "38")
    Output: list of dictionaries with keys: 
    "start", "end", "strand", "name", "description", "biotype", "ID", "distance", "orientation" 
    '''

    end=pos+window
    start=pos-window
    if start<1: 
        start=1

    overlapping_genes=restQuery(makeGeneOverlapQueryURL(chrom,start,end,build=build))
    gene_list=list()

    for gene in overlapping_genes:
        d_from_start=abs(pos-int(gene["start"]))
        d_from_end=abs(pos-int(gene["end"]))
        distance=min(d_from_start,d_from_end)
        if pos>=int(gene["start"]) and pos<=int(gene["end"]):
            distance=0

        gene_details = {
            "start" : int(gene["start"]),
            "end" : int(gene["end"]),
            "strand" : gene["strand"],
            "name" : gene["external_name"],
            "description" : gene["description"],
            "biotype" : gene["biotype"],
            "ID" : gene["id"],
            "distance":distance
        }

        # orientation:
        gene_details["orientation"] = "upsteram"
        if distance==0:
            gene_details["orientation"] = "overlapping"
        elif gene["strand"] == 1 and d_from_end < d_from_start:
            gene_details["orientation"] = "downstream"
        elif gene["strand"] == -1 and d_from_end > d_from_start:
            gene_details["orientation"] = "downstream"

        gene_list.append(gene_details)

    return gene_list

# ======================================================================================================================

# Extract Exome Aggregation Consortium (ExAC) allele frequencies
def getExacAF(chrom,pos,ref,alt):
    '''
    For a given variant, retreive allele frequencies for all ALT alleles in all populations
    Input  : chromosome, position, REF, ALT
    Output : list of dictionaries with the following keys "allele", "populations"
    "populations" value is a dictionary with "count" and "frequency" keys
    '''

    pops=["AFR", "ALL", "FIN", "AMR", "EAS", "SAS", "OTH", "NFE"]

    ExAC_file=config.EXAC_FILE
    if not os.path.isfile(ExAC_file):
        LOGGER.error("ExAC file %s not found" %(ExAC_file))
        return None

    query="tabix %s %s:%s-%s" %(ExAC_file,chrom,str(pos),str(pos))
    output=subprocess.Popen(query.strip(),shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    parsed=list()
    for line in output.stdout.readlines():
        fields=line.strip().split("\t")
        alts = fields[4].split(",")
        if (fields[0]==chrom and int(fields[1])==pos and fields[3]==ref and alt in alts):
            data=dict()
            for x in fields[7].split(";"):
                try:
                    (key, value)=x.split("=")
                    data[key]=value.split(",")
                except:
                    pass

            for i in range(0,len(alts)):
                d=dict()
                d["allele"]=alts[i]
                d["populations"]=dict()
                for p in pops:
                    label1="AC_"+p
                    label2="AN_"+p
                    if p=="ALL": 
                        label1="AC"
                        label2="AN"
                    count="NA"
                    freq="NA"
                    if label1 in data and label2 in data:
                        count=int(data[label1][i])
                        freq=float(count)/int(data[label2][0])
                    d["populations"][p]={"count":count,"frequency":freq}
                parsed.append(d)
    
    return parsed

# ======================================================================================================================

# Retrieve UK10K allele frequency:
# def get_UK10K_frequencies(variant_data):
#     '''
#     This function queries the uk10k public data for allele frequency and allele counts
#     It need the previously generate variant_data dictionary with all details:
#         - chromosome, variant location, start, end alternative, reference, allele string, rsID

#     Returned data: dictionary with the following structure:
#     {'Allele_counts': {
#         'A': 7446,
#         'G': 116},
#      'Allele_frequencies': {
#         'A': 0.98466014281936,
#         'G': 0.015339857180640043},
#      'Genotype_counts': {
#         'AA': '3665 (3665.89)',
#         'GA': '116 (114.22)',
#         'GG': '0 (0.89)'}}
#     In the genotype counts, the number in parentheses indicates the expected counts assuming HWE

#     '''

#     chromosome = variant_data["chromosome"]
#     queried_var = variant_data["location"]
#     start = variant_data["start"]
#     end = variant_data["end"]
#     alternative = variant_data['alternativeAllele']
#     reference  = variant_data["reference"]
#     allele_string = variant_data["allele_string"]
#     rsID = variant_data["end"]

#     # Check if datafile is exists, and return 0 if not:
#     return "UK10K datafiles were not found."

#     # tabix indexed vcf file, in which we are looking for the variation:
#     UK10K_vcf = config.UK10K_FILE % (chromosome)

#     # Check if datafile is exists, and return 0 if not:
#     if not os.path.isfile(UK10K_vcf):
#         return "UK10K datafiles were not found."

#     # submit bcftools query:
#     query = "bash -O extglob -c \'/software/hgi/pkglocal/bcftools-1.2/bin/bcftools query -f \"%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n\" -r %s %s\'" %(queried_var, UK10K_vcf)
#     output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

#     # Initializing returned variable:
#     genotype_count = {}
#     for line in output.stdout.readlines():

#         # from each lines, the info fields are extracted for checking:
#         fields = line.strip().split("\t")
#         out_chromosome = fields.pop(0)
#         out_position = fields.pop(0)
#         out_ID = fields.pop(0)
#         out_reference = fields.pop(0)
#         out_alternative = fields.pop(0)

#         # Now we test if the variation was found in the UK10K database:
#         if rsID == out_ID or (chromosome == out_chromosome and
#             (int(out_position) == start or int(out_position) == end) and
#             out_alternative in allele_string and out_reference in allele_string):

#             # Now we parse the genotypes:
#             genotype_count = {
#                 "Allele_counts": {
#                     out_reference : 0,
#                     out_alternative : 0
#                 },
#                 "Genotype_counts" : {
#                     out_reference+out_reference : 0,
#                     out_alternative+out_reference : 0,
#                     out_alternative+out_alternative : 0,
#                 },
#                 "Allele_frequencies": {
#                     out_reference : 0,
#                     out_alternative : 0
#                 }
#             }

#             n = 0 # Counting samples in row
#             for sample in fields:
#                 try:
#                     gt = sample.split("=")[1]
#                     (a1, a2) = gt.split("|")

#                     n += 1

#                     # Parsing genotypes:
#                     if a1 == "1" and a2 == "1":
#                         genotype_count["Allele_counts"][out_alternative] += 2
#                         genotype_count["Genotype_counts"][out_alternative+out_alternative] += 1
#                     if a1 == "0" and a2 == "0":
#                         genotype_count["Allele_counts"][out_reference] += 2
#                         genotype_count["Genotype_counts"][out_reference+out_reference] += 1
#                     else:
#                         genotype_count["Allele_counts"][out_reference] += 1
#                         genotype_count["Allele_counts"][out_alternative] += 1
#                         genotype_count["Genotype_counts"][out_alternative+out_reference] += 1

#                 except:
#                     pass # Unexpected field. Skip

#             # Calculating the allele frequencies:
#             genotype_count["Allele_frequencies"][out_reference] = round(float(genotype_count["Allele_counts"][out_reference]) /(
#                                     genotype_count["Allele_counts"][out_alternative] + genotype_count["Allele_counts"][out_reference]), 3)
#             genotype_count["Allele_frequencies"][out_alternative] = round(float(genotype_count["Allele_counts"][out_alternative]) / (
#                                     genotype_count["Allele_counts"][out_alternative] + genotype_count["Allele_counts"][out_reference]), 3)

#             # Calculating expected number of genotypes (based on pop size and HW):
#             exp_AA = (genotype_count["Allele_frequencies"][out_reference] ** 2) * n
#             exp_aa = (genotype_count["Allele_frequencies"][out_alternative] ** 2) * n
#             exp_Aa = 2 * (genotype_count["Allele_frequencies"][out_alternative] * genotype_count["Allele_frequencies"][out_reference]) * n

#             # Updating genotype counts with the expected values:
#             genotype_count["Genotype_counts"][out_alternative+out_alternative] = '%s (%s)'  %(genotype_count["Genotype_counts"][out_alternative+out_alternative],
#                                                                                            round(exp_aa, 2))
#             genotype_count["Genotype_counts"][out_reference+out_reference] = '%s (%s)'  %(genotype_count["Genotype_counts"][out_reference+out_reference],
#                                                                                           round(exp_AA, 2))
#             genotype_count["Genotype_counts"][out_alternative+out_reference]= '%s (%s)'  %( genotype_count["Genotype_counts"][out_alternative+out_reference],
#                                                                                           round(exp_Aa, 2))

#             # If we found only one matching line, we return data:
#             return genotype_count

#         else:
#             pass # Going to the next line.


#     # Once we reached the end of the lines without finding the correct line, then we
#     return "Variation was not found in the UK10K databaset."

# ======================================================================================================================

def getRegulation(chrom,pos,window=2000):
    '''
    This function returns all overlapping active regulatory features within window of the variation.

    Input: chromosome, position, window (default: 2000)

    Output: dictionary with all regulatory IDs as keys, dictionary("chrom","start","end","class","cells") as values

    '''

    # Cell type abbreviations are retrieved from Ensembl, which are replaced
    # by conventional names for better readibility.
    CellTypeFinder = config.CellTypeFinder

    start=pos-window
    if start<1:
        start=1
    end=pos+window

    regulatoryFile = config.REGULATORY_FILE
    FNULL = open(os.devnull, 'w')
    query="intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\n\") -b %s -sorted" % (chrom,start,end,regulatoryFile)
    output=subprocess.Popen(query.strip(),shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=FNULL)

    parsed=dict()
    for line in output.stdout.readlines():
        #print(line)
        c0,s0,e0,c1,st,en,cl,cell,act,regid=line.strip().split("\t")
        if act=="ACTIVE":
            if regid not in parsed:
                parsed[regid]={"chrom":c1,"start":int(st),"end":int(en),"class":cl,"cells":[cell]}
            else:
                if parsed[regid]["chrom"]!=c1 or parsed[regid]["start"]!=int(st) or parsed[regid]["end"]!=int(en) or parsed[regid]["class"]!=cl:
                    LOGGER.error("Regulatory ID %s has conflicting attributes: %s %d %d %s (%s %d %d %s)" % (regid,parsed[regid]["chrom"],parsed[regid]["start"],parsed[regid]["end"],parsed[regid]["class"],c1,int(st),int(en),cell))
                    FNULL.close()
                    return None
                else:
                    if cell not in parsed[regid]["cells"]:
                        parsed[regid]["cells"].append(cell)

    FNULL.close()
    return parsed

# ===================================== CONVERTING DIFFERENT DATA STRUCTURES TO DATAFRAMES =====================================

# TODO: add links like in the original version
# TODO: gwava and gerp
def variant2df(var_data):
    df=pd.DataFrame(columns=["Value"])
    df.loc["ID"]=[var_data["rsID"]]
    for m in var_data["mappings"]:
        df.loc["Location"]=[getLocationString(m["chr"],m["pos"],m["ref"],m["alt"])]
        df.loc["Allele string"]=[m["ref"]+"/"+m["alt"]]
        df.loc["SIFT score"]=[m["sift_score"]]
        df.loc["SIFT prediction"]=[m["sift_prediction"]]
        df.loc["PolyPhen score"]=[m["polyphen_score"]]
        df.loc["PolyPhen prediction"]=[m["polyphen_prediction"]]
    df.loc["MAF"]=[var_data["MAF"]]
    df.loc["Consequence"]=[var_data["consequence"]]
    df.loc["Type"]=[var_data["class"]]
    s=", ".join(var_data["synonyms"])
    if len(var_data["synonyms"])>5:
        s=",".join(var_data["synonyms"][0:5])
        s=s+"... ("+str(len(var_data["synonyms"]))+" in total)"
    df.loc["Synonyms"]=[s]
    for p in var_data["phenotype_data"]:
        df.loc["Phenotype"]=[p["trait"]+" - "+p["source"]+" - "+p["risk_allele"]]
    if len(var_data["clinical_significance"])!=0:
        df.loc["Clinical significance"]=[",".join(var_data["clinical_significance"])]

    return df
    
# ----------------------------------------------------------------------------------------------------------------------

def regulation2df(reg_data):
    df=pd.DataFrame(columns=["ID","Class","Cell type"])
    i=0
    for r in reg_data:
        for cell in reg_data[r]["cells"]:
            df.loc[i]=[r,reg_data[r]["class"],cell]
            i+=1
    
    return df

# ----------------------------------------------------------------------------------------------------------------------

def gwas2df(gwas_data):
    df=pd.DataFrame(columns=["rsID","ID","Distance","Trait","P-value","PMID"])
    i=0
    for x in gwas_data:
        df.loc[i]=[x["rsID"],x["SNPID"],x["distance"],x["trait"],x["p-value"],x["PMID"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def vepTranscript2df(vep_data):
    df=pd.DataFrame(columns=["Gene ID","Transcript ID","Impact","Consequence","Isoform"])
    i=0
    for x in vep_data["transcript"]:
        df.loc[i]=[x["gene_id"],x["ID"],x["impact"],x["consequence"],x["principal"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def vepRegulatory2df(vep_data):
    df=pd.DataFrame(columns=["Biotype","Regulatory feature ID","Impact","Consequence"])
    i=0
    for x in vep_data["regulatory"]:
        df.loc[i]=[x["biotype"],x["ID"],x["impact"],",".join(x["consequence"])]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def population2df(pop_data):
    df=pd.DataFrame(columns=["Population","Allele 1","Allele 2"])
    i=0
    for p in config.PopulationNames:
        x=next((z for z in pop_data if z["population"]==p),None)
        if x:
            L=list(x["frequency"])
            L.sort()
            df.loc[i]=[p,L[0]+" ("+str(round(x["frequency"][L[0]],4))+")",L[1]+" ("+str(round(x["frequency"][L[1]],4))+")"]
            i+=1
        else:
            df.loc[i]=[p,"NA","NA"]
            i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def pubmed2df(pubmed_data):
    df=pd.DataFrame(columns=["First author","Journal","Year","URL","Title"])
    i=0
    for x in pubmed_data:
        d=pubmed_data[x]
        df.loc[i]=[d["firstAuthor"],d["journal"],d["year"],"<a href='"+d["URL"]+"'>Link</a>",d["title"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def geneList2df(gene_data):
    df=pd.DataFrame(columns=["Gene name","Gene ID","Biotype","Distance","Orientation"])
    i=0
    for x in gene_data:
        df.loc[i]=[x["name"],x["ID"],x["biotype"],x["distance"],x["orientation"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def exac2df(exac_data):
    L=["Population"]
    L1=[]
    for x in exac_data:
        L1.append(x["allele"])
    L1.sort()
    for a in L1:
        L.append(a)

    df=pd.DataFrame(columns=L)
    ExAC_pops = ['NFE','FIN','AFR','AMR','EAS','SAS','OTH','ALL']
    i=0
    for p in ExAC_pops:
        L=[p]
        for a in L1:
            x=next((z for z in exac_data if z["allele"]==a),None)
            if x:
                if p in x["populations"]:
                    c=x["populations"][p]["count"]
                    f=x["populations"][p]["frequency"]
                    L.append("%s (%s)" %(str(c),str(round(f,4))))
                else:
                    L.append("NA (NA)")
            else:
                LOGGER.error("Could not find allele %s" %(a))

        df.loc[i]=L
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def gtex2df(gtex_data):
    df=pd.DataFrame(columns=["Tissue","P-value","Beta (SE)","ID","Distance from TSS"])
    i=0
    for x in gtex_data:
        for z in gtex_data[x]:
            df.loc[i]=[z["tissue"],z["p-value"],z["beta"]+" ("+z["SE"]+")",x,z["dist"]]
            i+=1

    return df

# =================================== CONVERTING GENE RELATED DATA STRUCTURES TO DATAFRAMES =============================

def geneInfo2df(gene_info):
    df=pd.DataFrame(columns=["Gene name","Description","ID","Coordinates","Strand","Type"])
    df.loc[0]=[gene_info["name"],gene_info["description"],gene_info["id"],gene_info["chromosome"]+":"+str(gene_info["start"])+"-"+str(gene_info["end"]),gene_info["strand"],gene_info["type"],]

    return df

# ----------------------------------------------------------------------------------------------------------------------

def uniprot2df(data):
    df=pd.DataFrame(columns=["Field","Data"])
    i=0
    for k in data:
        df.loc[i]=[k,data[k]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def geneGwas2df(data):
    df=pd.DataFrame(columns=["rsID","Consequence","Reported genes","Trait","URL","Author"])
    i=0
    for x in data:
        df.loc[i]=[x["rsID"],x["Consequence"],x["Reported genes"],x["trait"],"<a href='https://"+x["URL"]+"'>Link</a>",x["Author"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def goterms2df(xrefs):
    df=pd.DataFrame(columns=["GO term ID","Description"])
    i=0
    for x in xrefs["GO"]:
        df.loc[i]=[x[0],x[1]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def generateHTML(templateFile, data):
    templateLoader = jinja2.FileSystemLoader(searchpath="/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(templateFile)
    return template.render(data)

