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

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

def makeOverlapVarQueryURL(chrom,start,end,build="38"):
    ext = "/overlap/region/human/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s:%s-%s?feature=variation"  % (chrom,str(start),str(end))

def makeRefSeqQueryURL(chrom,start,end,build="38"):
    ext = "/sequence/region/human/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s:%s..%s:1?"  % (chrom,str(start),str(end))

def makeRSQueryURL(rsID,build="38"):
    ext = "/variant_recoder/homo_sapiens/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+rsID+"?"

def makeHomologyURL(ID,species="mouse",build="38"):
    ext="/homology/id/%s?&target_species=%s&aligned=0&sequence=none&type=orthologues" %(ID,species)
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext

def makeGeneQueryURL(ID,build="38"):
    ext="/lookup/id/"
    server="http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+ID

def makeGeneXQueryURL(ID,build="38"):
    ext="/xrefs/id/"
    server="http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s?all_levels=1" % ID

def makeGeneXQueryURL2(ID,build="38"):
    ext="/xrefs/id/"
    server="http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s?external_db=MGI" % ID

def makeRSPhenotypeQueryURL(build="38"):
    ext = "/variation/homo_sapiens/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"?phenotypes=1"

def makePhenoOverlapQueryURL(chrom,start,end,build="38"):
    ext = "/overlap/region/human/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s:%d-%d?feature=variation;variant_set=ph_variants;content-type=application/json" %(chrom,start,end)

def makeRsPhenotypeQuery2URL(rs,build="38"):
    ext = "/variation/human/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+"%s?pops=1;phenotypes=1" %rs

def makeRSListQueryURL(build="38"):
    ext = "/variant_recoder/homo_sapiens/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext

def makeVepQueryURL(chrom,start,end,allele,strand="1",build="38"):
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"
    ext="/vep/human/region/%s:%s-%s:1/%s?" % (chrom,str(start),str(end),allele)

    return server+ext

#------------------------------------------------------------------------------------------------------------------------------------

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
                    print("[WARNING]: parseSPDI: deletion: (%s %d %s %s) INS is not suffix of DEL" %(c,pos,ref,alt),file=sys.stderr)
                    return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
            # insertion
            elif len(delseq)<len(alt):
                if lref==0:
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt}
                elif alt.endswith(delseq):
                    alt2=alt[0:len(alt)-len(ref)]
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt2}
                else:
                    print("[WARNING]: parseSPDI: insertion: (%s %d %s %s) DEL is not prefix of INS" %(c,pos,ref,alt),file=sys.stderr)
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
        print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: query type ("+qtype+") has to be either \"get\" or \"post\"",file=sys.stderr)
        sys.stderr.flush()
        return None

    r=None
    try:
        if qtype=="get":
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},timeout=timeout)
        else:
            if not data:
                print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: Error: postquery requires data",file=sys.stderr)
                sys.stderr.flush()
                return None                
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},data=data,timeout=timeout)

        if not r.ok:
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: Error "+str(r.status_code)+" occured (input: "+URL+")",file=sys.stderr)
            sys.stderr.flush()
            return None

        try:
            ret=r.json()
            return ret
        except ValueError:
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: JSON decoding error", file=sys.stderr)
            sys.stderr.flush()
            return None

    except Timeout as ex:
        print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: Timeout exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except TooManyRedirects as ex:
        print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: TooManyRedirects exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except RequestException as ex:
        print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: RequestException occured", file=sys.stderr)
        sys.stderr.flush()
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


# for a given genomic region, return dataframe containing variants with phenotype annotations
def getVariantsWithPhenotypes(chrom,start,end,pos=0,build="38"):
# 5Mb max !

    '''
    output: pandas dataframe with the following columns:
    'SNPID', 'consequence', 'distance', 'genes', 'phenotype', 'rsID', 'source'
    '''
    if start<0 : start=0
    URL=makePhenoOverlapQueryURL(chrom,start,end,build=build)

    variants=restQuery(URL,qtype="get")

    if not variants:
        return None

    if len(variants) == 0: 
        print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getVariantsWithPhenotypes: No variants with phenotypes were found in the region", file=sys.stderr)
        return None

    rsIDs = []
    for var in variants:
        rsIDs.append(var["id"])

    # Given rsIDs, get the phenotypes
    URL=makeRSPhenotypeQueryURL(build=build)
    r = restQuery(URL,data=list2string(rsIDs),qtype="post")

    # Check output:
    if not r:
        return None

    buildstr=""
    if build != "38":
        buildstr="grch"+build+"."

    output=[]
    for rsID in r:
        for phenotype in r[rsID]["phenotypes"]:

            # TODO: check all possible sources
            # Generate phenotype link based on source:
            if phenotype["source"] == "ClinVar":
                link = "https://www.ncbi.nlm.nih.gov/clinvar/?term="+rsID
            elif phenotype["source"] == "HGMD-PUBLIC":
                link = "http://www.hgmd.cf.ac.uk/ac/gene.php?gene="+phenotype["genes"]
            elif "NHGRI-EBI" in phenotype["source"]:
                link = "https://www.ebi.ac.uk/gwas/search?query="+rsID
            elif phenotype["source"] == "OMIM":
                link = "https://www.omim.org/entry/"+phenotype["study"][4:]
            else:
                link = "http://"+buildstr+"ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;v=CM106680;vdb=variation"

            # TODO: check key availability
            output.append({"rsID": rsID,
                    "consequence" : r[rsID]["most_severe_consequence"],
                    "SNPID" : "%s:%s" %(
                        r[rsID]["mappings"][0]["seq_region_name"],
                        r[rsID]["mappings"][0]["start"]),
                    "phenotype" : phenotype["trait"],
                    "URL" : link,
                    "distance" : abs(pos - r[rsID]["mappings"][0]["start"])})

    # Create a dataframe:
    df = pd.DataFrame(output)
    df.consequence = df.consequence.str.title().str.replace("_", " ")

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
    res["class"]=data["var_class"]
    res["consequence"]=data["most_severe_consequence"]

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
            mappings.append({"chr":c,"pos":p,"ref":ref,"alt":alt})

#------------------ population data ----------------

    population_data=list()

    for pop in data["populations"]:
        pop_name = pop["population"].split(":")
        if pop_name[0] == "1000GENOMES" and pop_name[1] == "phase_3":
            name=pop_name[2]
            try:
                z=next(x for x in population_data if name in x.values())
                z[pop["allele"]]=pop["frequency"]
            except:
                population_data.append({"population":name,pop["allele"]:pop["frequency"]})

#------------------ phenotype data -------------------

    phenotype_data=list()

    for p in data["phenotypes"]:
        trait=p["trait"] if "trait" in p else None
        source=p["source"] if "source" in p else None
        risk=p["risk_allele"] if "risk_allele" in p else None
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

    m=re.match("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
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
        for v in r:
            z=restQuery(makeRSQueryURL(v["id"],build=build))
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

                    if (ref==a1 and alt==a2) or (ref==a2 and alt==a1)
                    #if seq0 in stringdiff(alt[1:len(alt)],ref[1:len(ref)]):
                        S.add(var)
                        break
    return S

# ===================================================== GTEx RELATED STUFF ============================================

# given gene ID (variant ID), retreive all variant (gene) data associated with the gene (variant): tissue, p-value
def parseGTEx(filename,chrom,start,end,ID):
    query = "bash -O extglob -c \'tabix %s %s:%d-%d | fgrep -w %s\'" %(filename,chrom,start,end,ID)
    output = subprocess.Popen(query.strip(),universal_newlines=True,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    D=dict()
    for line in output.stdout.readlines():
        fields=line.strip().split("\t")
        data=fields[4].split(":")
        ID2=data[0]
        tissue=data[1]
        pval=data[2]
        if not ID2 in D:
            D[ID2]=[]
        D[ID2].append({"tissue" : tissue,"pvalue" : pval})

    if len(D)==0:
        print("INFO: no eQTL signals were found for %s" %(ID),file=sys.stderr)
    
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

# retrieve a list of gwas signals around the variant
def getGwasHits(chrom,pos,window=500000):
    '''
    This function retrives a list of all gwas signals within 500kbp distance around a given position.

    input: chromosome, position

    output: [ {"rsID","SNPID","trait","p-value","PMID","distance"}]
    '''

    gwas_file = config.GWAS_FILE_VAR

    if not os.path.isfile(gwas_file):
        print("[ERROR] GWAS catalog file (%s) not found" % gwas_file,file=sys.stderr)
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
    This function retrieve gene related information

    INPUT: Ensembl stable ID
    OUTPUT: dictionary with retrieved information
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
    '''
    response=restQuery(makeGeneXQueryURL(ID,build=build))
    xrefs = {
        "MIM disease" : [],
        "MIM gene" : [],
        "GO" : [],
        "GOSlim GOA" : [],
        "UniProtKB/Swiss-Prot" : [],
        "Human Protein Atlas" : [],
        "ChEMBL" : [],
    }

    for xref in response:
        db=xref['db_display_name']
        try:
            xrefs[db].append([xref["description"],xref["primary_id"]])
        except:
            continue

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
        print("[INFO] No mouse cross-reference for %s" %(human_ID),file=sys.stderr)
        return mouse_IDs

    for homolog in data["data"][0]["homologies"]:
        try:
            z=restQuery(makeGeneQueryURL(homolog["target"]["id"],build=build))
            #print(json.dumps(z,indent=4,sort_keys=True))

            d=""
            m=re.match("(.*)\s+\[.*\]",z["description"])
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

    print("[INFO] MGI ID for %s was not found" % (mouse_ID),file=sys.stderr)
    return ""

def getMgiPhenotypes(MGI_ID):
    '''Returns phenotype information stored on http://www.informatics.jax.org/ '''

    URL = "http://www.informatics.jax.org/allele/report.txt?markerId=%s" % MGI_ID

    # Return all associated allele information:
    try:
        df=pd.read_csv(URL, sep="\t")

        # Drop unnecessary columns:
        df.drop([ "Allele Symbol", "Chromosome",  "Synonyms","Allele Attributes", "Transmission", "Unnamed: 10"], axis=1, inplace=True)
        df.columns = ["Allele_ID", "Allele_name", "Allele_type","Phenotypes","Human_disease"]
        return df
    except:
        print("[INFO] no phenotype was found for %s" %(MGI_ID),file=sys.stderr)
        return pd.DataFrame(columns=["Allele_ID", "Allele_name", "Allele_type","Phenotypes","Human_disease"])

def getMousePhenotypes(geneID):
    '''returning mouse phenotype given human gene ID'''
    # Returning all mouse homologue IDs:
    mouse_gene_IDs=getMouseID(geneID)

    full_dataframe=pd.DataFrame(columns=["Allele_ID", "Allele_name", "Allele_type","Phenotypes","Human_disease","mouse_gene_ID","MGI_ID","mouse_gene_name","mouse_gene_description"])

    if len(mouse_gene_IDs)==0:
        return full_dataframe

    MGI_IDs=dict()
    for mouse_gene_ID in mouse_gene_IDs:
        MGI_IDs[mouse_gene_ID] = getMgiID(mouse_gene_ID)

    # Once we have all the MGI identifiers, we retrieve all the phenotypes:
    for mouse_id, mgi_id in MGI_IDs.items():
        df=getMgiPhenotypes(mgi_id)

        # Adding extra columns for the record:
        df["mouse_gene_ID"] = mouse_id
        df["MGI_ID"] = mgi_id
        df["mouse_gene_name"] = mouse_gene_IDs[mouse_id][1]
        df["mouse_gene_description"] = mouse_gene_IDs[mouse_id][0]
        full_dataframe = pd.merge(full_dataframe, df, how="outer")

    return full_dataframe

# ======================================================= APPRIS ===========================================================

def getApprisInfo(gene_ID):
    '''
    This function returns transcript information from the APPRIS webserver based on Ensembl gene ID
    The output is a dictionary with detailed information of all annotated transcripts

    Link: http://appris.bioinfo.cnio.es/
    '''

    # hardcoded assembly
    URL = "http://apprisws.bioinfo.cnio.es:80/rest/exporter/id/homo_sapiens/%s?format=json&db=b38" % gene_ID;

    data=dict()
    r = requests.get(URL)
    if not r.ok:
        print("[INFO] Query failed for gene: %s\nURL: %s" % (gene_ID,URL),file=sys.stderr)
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

# Retiveing data from the variant effect predictor server.
def getVepData(variant_data):
    '''
    This function returns the predicted effect based on chromosome, position and the alternate allele.

    Input: variant data, a dictionary with chr, pos, ref, alt

    Output: variant effects
    '''

    chrom = variant_data["chr"]
    pos = variant_data["pos"]
    ref = variant_data["ref"]
    alt = variant_data["alt"]

    start=-1
    end=-1
    allele=""

    if len(ref)>1 and len(alt)==1:
        start=pos+1
        end=pos+len(ref)-1
        allele=ref[1:]
    elif len(alt)>1 and len(ref)==1:
        start=pos+1
        end=pos
        allele=alt[1:]
    elif len(ref)==1 and len(alt)==1:
        start=pos
        end=pos
        allele=alt
    else:
        print("[ERROR]: wrong allele encoding ref=%s, alt=%s" %(ref,alt),file=sys.stderr)

    data=restQuery(makeVepQueryURL(chrom,start,end,allele))
    return data

    # The most severe Sift and Polyphen scores will be retrieved:
    sift = []
    polyphen = []

    # Big dataset for all transcripts and regulatory features will be retrieved:
    VEP_data = {'transcript' : [],
                'regulatory' : []}

    # We retrieve the amino acid changes as well:
    aminoacid = ""
    codon = ""

    # If there is no consequence in the variant data then, we add: <-It does not work....
    variant_data["consequence"] = VEP[0]['most_severe_consequence']

    # Gathering transcript consequences:
    if 'transcript_consequences' in VEP[0].keys():
        VEP_data['transcript']=[]
        for transcript in VEP[0]['transcript_consequences']:
            VEP_data['transcript'].append({
                "impact" : transcript["impact"],
                "gene_symbol" : transcript["gene_symbol"],
                "gene_id": transcript["gene_id"],
                "transcript_id": transcript["transcript_id"],
                "consequence" : transcript["consequence_terms"]
            })

            # Testing polyphen score:
            if "polyphen_score" in transcript.keys():
                try:
                    if transcript["polyphen_score"] > polyphen[0]:
                        polyphen[0] = transcript["polyphen_score"]
                        polyphen[1] = transcript["polyphen_prediction"]
                except:
                        polyphen.append(transcript["polyphen_score"])
                        polyphen.append(transcript["polyphen_prediction"])

            # Testing sift score:
            if "sift_score" in transcript.keys():
                try:
                    if transcript["sift_score"] < sift[0]:
                        sift[0] = transcript["sift_score"]
                        sift[1] = transcript["sift_prediction"]
                except:
                        sift.append(transcript["sift_score"])
                        sift.append(transcript["sift_prediction"])

            # Testing amino acid substitution:
            if VEP[0]["most_severe_consequence"] in transcript["consequence_terms"] and "amino_acids" in transcript.keys():
                aminoacid = transcript["amino_acids"].replace("/", str(transcript["protein_start"]))
                codon = transcript["codons"]


    # Gathering regulatory consequences:
    if 'regulatory_feature_consequences' in VEP[0].keys():
        VEP_data['regulatory']=[]
        for regulatory in VEP[0]['regulatory_feature_consequences']:
            VEP_data['regulatory'].append({
                "impact" : regulatory["impact"],
                "biotype" : regulatory["biotype"],
                "gene_id": "-",
                "regulatory_ID": regulatory["regulatory_feature_id"],
                "consequence" : regulatory["consequence_terms"]
            })

    # Updating variation data:
    variant_data["sift"] = sift
    variant_data["polyphen"] = polyphen
    variant_data["codon"] = codon
    variant_data["aminoacid"] = aminoacid

    # If there were overlapping transcripts, we have to check if the
    # get a list of all overlapping genes:
    gene_list = []
    for transcripts in VEP_data['transcript']:
        if not transcripts["gene_id"] in gene_list:
            gene_list.append(transcripts["gene_id"])

    # Now looping trugh all genes and get list of principal transcripts:
    for gene_id in gene_list:
        appris_data = get_principal_transcript(gene_id) # downloading appris data

        # Now looping trough all consequences to get the appris annotation:
        for transcripts in VEP_data['transcript']:
            if transcripts["gene_id"] in gene_id:
                try:
                    transcripts["principal"] = appris_data[transcripts["transcript_id"]]["reliability"]
                except:
                    transcripts["principal"] = "NA"

    # Checking if there is anything in the VEP data:
    if len(VEP_data['transcript']) == 0:
        VEP_data['transcript'] = "The queried variation does not overlap any transcripts."
    if len(VEP_data['regulatory']) == 0:
        VEP_data['regulatory'] = "The queried variation does not overlap any regulatory features."

    # Returning data
    return (VEP_data, variant_data)
