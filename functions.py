import requests
import datetime
import sys
import json 
import re
from requests.exceptions import Timeout,TooManyRedirects,RequestException

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

# This function downloads variation data from ensembl:
def get_variant_info(rsID, ref):
    '''
    This function returns the most important details of a variation given its rsID.

    It also returns the 1000 Genome Phase 3 allele frequencies as well.

    Returned: pop_freq, variation_data
    '''

    URL = config.REST_URL + "/variation/human/%s?content-type=application/json;pops=1;phenotypes=1" %rsID

    variationData = submit_REST(URL)

    # If the user only provides rsID, we don't know which is the reference allele:
    if ref == "NA":
        ref = variationData["ancestral_allele"]

    # Building our custom hash:
    variation_data = {
        "reference": ref,
        "rsID" : variationData["name"],
        "minorAllele" : variationData["minor_allele"],
        "ancestralAllele" : variationData["ancestral_allele"],
        "synonyms" : variationData["synonyms"],
        "MAF" : variationData["MAF"],
        "varClass" : variationData["var_class"],
        "consequence" : variationData["most_severe_consequence"]
    }

    # Extracting frequencies:
    pop_freq = {}
    for pops in variationData["populations"]:
        pop_name = pops["population"].split(":")
        try:
            if pop_name[0] == "1000GENOMES" and pop_name[1] == "phase_3":
                try:
                    pop_freq[pop_name[2]][pops["allele"]] = pops["frequency"]
                except:
                    pop_freq[pop_name[2]] = {pops["allele"] : pops["frequency"]}
        except:
            continue
    # If not 1000 genome frequencies are available for the given variation:
    if len(pop_freq) == 0:
        pop_freq = "Variation was not found in the 1000 Genomes data."

    # Extracting mapping:
    for mapping in variationData["mappings"]:
        if len(mapping["seq_region_name"]) > 2:
            continue
        else:
            variation_data["chromosome"] = mapping["seq_region_name"]
            variation_data["start"] = mapping["start"]
            variation_data["end"] = mapping["end"]
            variation_data["allele_string"] = mapping["allele_string"]
            variation_data["location"] = mapping["location"]

    # 'alternativeAllele': alt,
    if variationData["minor_allele"] == ref:
        all_alleles = variation_data["allele_string"].split("/")
        variation_data["alternativeAllele"] = all_alleles[1]
    else:
        variation_data["alternativeAllele"] = variationData["minor_allele"]

    # Extracting phenotype information:
    if len(variationData["phenotypes"]) > 0: variation_data["phenotypes"] = []
    for pt in variationData["phenotypes"]:
        trait = pt["trait"] if "trait" in pt else "-"
        source = pt["source"] if "source" in pt else "-"
        allele = pt["risk_allele"] if "risk_allele" in pt else "-"

        if "risk_allele" in pt:
            trait += "("+pt["risk_allele"]+")"

        variation_data["phenotypes"].append({
            "trait": trait,
            "source": source,
            })



    # Generating SNP ID:
    variation_data["SNP_ID"] = "chr" + str(variation_data["chromosome"]) + ":" + str(variation_data["start"])

    # Extracting clinical significance:
    try:
        for cs in variationData["clinical_significance"]:
            if cs != "other" and cs != "not provided":
                try:
                    variation_data["clinical_significance"].append(cs)
                except:
                    variation_data["clinical_significance"] = [cs]
    except:
        pass


    # Now we have to generate the gerp and gwava scores.
    # variation_data = get_GWAVA_score(variation_data)

    return (pop_freq, variation_data)

# return a part of the reference genome
def getRefSeq(chrom,start,end,build="38"):
    r=restQuery(makeRefSeqQueryURL(chrom,start,end,build))
    if r:
        return response["seq"]
    else:
        return None

# given variant ID, try to get variant type
# also check if reference sequence and alleles are not in contradiction
def getVariantType(varid,build="38"):
    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
    chrom=m.group(1)
    pos=int(m.group(2))
    a1=m.group(3)
    a2=m.group(4)

    vartype=None
    # check if reference sequence and alleles are not in contradiction
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

# returns a list of matching rsIDs for a given variant ID
def id2rs(varid,build="38"):
    if varid.startswith("rs"):
        return varid

    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
    chrom=m.group(1)
    pos=int(m.group(2))
    a1=m.group(3)
    a2=m.group(4)

    vartype=getVariantType(varid,build=build)

    # in case of indels, pull all variants around the given variant
    window=5

    L=[]
    if len(a1)==1 and len(a2)==1: # SNPs
        r=restQuery(makeOverlapVarQueryURL(chrom,pos,pos,build=build))
        for v in r:
            if a1 in v["alleles"] and a2 in v["alleles"]:
                L.append(v["id"])
    else:
        # difference between a1 and a2
        seq0=""
        if len(a1)>len(a2):
            seq0=a1[len(a2),len(a1)]
        elif len(a1)<len(a2):
            seq0=a2[len(a1),len(a2)]

        r=restQuery(makeOverlapVarQueryURL(chrom,pos-window,pos+window,build=build))
        for v in r:
            z=restQuery(makeRSQueryURL(v["id"],build=build))
            for x in z:
                spdi=x["spdi"]
                var=x["id"]
                h=parseSPDI(spdi,alleles=True)
                ref=h["ref"]
                alt=h["alt"]
                pos=h["pos"]
                if len(ref)>len(alt): # deletion
                    seq=""
                    if ref.startswith(alt):
                        seq=ref[len(alt),len(ref)] # ref's prefix
    return L
