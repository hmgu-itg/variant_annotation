import requests
import datetime
import sys
import json 
import re
import subprocess
from requests.exceptions import Timeout,TooManyRedirects,RequestException
import config

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
    m=re.search("^(\d+)$",ref)
    if m:
        isNum=True
        if ref=="1" and lalt==1: # one base deleted, one inserted, i.e. variant is a SNP
            isSNP=True
    elif lref==1 and lalt==1:# SNP
        isSNP=True

    if alleles:
        base=getRefSeq(c,pos,pos,build)

        # deleted sequence
        if isNum:
            seq=getRefSeq(c,pos+1,pos+int(ref),build)
        else:
            seq=ref

        if isSNP:
            return {"chr":c,"pos":pos+1,"ref":seq[0],"alt":alt}
        else:
            return {"chr":c,"pos":pos,"ref":base+seq,"alt":base+alt}
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
            print(str(datetime.datetime.now().strftime("%H:%M:%S"))+" : getQuery: Error "+str(r.status_code)+" occured",file=sys.stderr)
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

# get general information about a variant, given rsID:
# 
# MAF, minor allele, variant class, most severe consequence
# 1KG phase 3 population allele frequencies
# mappings: chr, pos, ref, alt
# phenotype data: trait, source, risk allele
# variant's clinical significance
def getVariantInfo(rs,build="b38"):
    res=dict()

#------------------- general information ---------------
    data=restQuery(makeRsPhenotypeQuery2URL(rs,build))

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

    for cs in data["clinical_significance"]:
        if cs != "other" and cs != "not provided":
            clinical_significance.append(cs)

#-----------------------------------------------------

    res["mappings"]=mappings
    res["population_data"]=population_data
    res["phenotype_data"]=phenotype_data
    res["clinical_significance"]=clinical_significance

    return res

# return a part of the reference genome
def getRefSeq(chrom,start,end,build="38"):
    r=restQuery(makeRefSeqQueryURL(chrom,start,end,build))
    if r:
        return r["seq"]
    else:
        return None

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

# returns a list of matching rsIDs for a given variant ID
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
    if len(a1)==1 and len(a2)==1: # SNPs
        r=restQuery(makeOverlapVarQueryURL(chrom,pos,pos,build=build))
        for v in r:
            if a1 in v["alleles"] and a2 in v["alleles"]:
                S.add(v["id"])
    else:
        # difference between a1 and a2
        seq0=""
        if len(a1)>len(a2):
            seq0=a1[len(a2):len(a1)]
        elif len(a1)<len(a2):
            seq0=a2[len(a1):len(a2)]

        r=restQuery(makeOverlapVarQueryURL(chrom,pos-window,pos+window,build=build))
        for v in r:
            z=restQuery(makeRSQueryURL(v["id"],build=build))
            for x in z:
                spdis=x["spdi"]
                var=x["id"][0]
                #print("")
                #print(varid)
                #print(json.dumps(x, indent=4, sort_keys=True))
                for spdi in spdis:
                    #print("SPDI: "+spdi)
                    h=parseSPDI(spdi,alleles=True)
                    ref=h["ref"]
                    alt=h["alt"]
                    p=h["pos"]
                    c=h["chr"]

                    #print("chr: %s, pos: %d, ref: %s, alt: %s" % (c,p,ref,alt))

                    if p!=pos:
                        continue

                    if len(ref)==1 and len(alt)==1:
                        continue
                        
                    if seq0 in stringdiff(alt[1:len(alt)],ref[1:len(ref)]):
                        S.add(var)
                        break
    return S

#===================================================== GTEx RELATED STUFF ============================================

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

#===================================================== UNIPROT RELATED STUFF ============================================

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

    r = requests.get(URL)
    print(r.text)
    print(r.encoding)

    fields=r.content.decode("utf-8").split("\n")[1].split("\t")

    UniprotData = {
        "Name" : fields[9],
        "Disease" : fields[8],
        "Function": fields[1],
        "Entry" : fields[0],
        "Subunit" : fields[2],
        "Phenotype" : fields[6],
        "localization" : fields[7],
        "Tissue" : fields[4],
        "Development" : fields[3]
    }

    return UniprotData


# Retrieve uniprot data:
def get_UNIPROT_data(cross_refs):
    '''
    The cross references are taken as input. It might contain none, one or
    multiple references to the uniprot database. The function loops through all
    and returns a dictionary with all unique protein names as keys.

    Retrieved fields:
        Function, Subunit, Developmental stage,
        Tissue specificity, Catalytic activity,
        Disruption phenotype, Subcellular localiztion,
        Disease

    The list can be extended if need....
    '''

    # Looping through all cross references pointing to Uniprot:
    Independent_IDs = {}
    for uniprotID in cross_refs["UniProtKB/Swiss-Prot"]:
        Independent_IDs[uniprotID[0]] = uniprotID[1]

    # Annotate all unique proteins:
    annotated_protein  = {}
    for ID in Independent_IDs.keys():
        annotated_protein[ID] = __UNIPROT_DOWNLOAD__(Independent_IDs[ID])

    # Returning annotation:
    return (annotated_protein)
