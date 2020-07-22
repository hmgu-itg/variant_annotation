import config
import logging
import pandas as pd
import json
import re
import os
import tempfile as tf
import subprocess

from query import *
from utils import *

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==============================================================================================================================

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

# for a given genomic region, return dataframe containing variants with phenotype annotations
def getVariantsWithPhenotypes(chrom,pos,window=config.WINDOW,build="38"):
    '''
    For a given genomic region, return dataframe containing variants with phenotype annotations

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

    LOGGER.debug("%s:%d" %(chrom,pos))
    variants=restQuery(makePhenoOverlapQueryURL(chrom,start,end,build=build),qtype="get")
    #print(json.dumps(variants,indent=4,sort_keys=True))

    if not variants:
        return None

    if len(variants)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return None

    rsIDs=list()
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
            #print(json.dumps(r,indent=4,sort_keys=True))
            for rsID in r:
                for phenotype in r[rsID]["phenotypes"]:
                    m=re.match(".*phenotype\s+not\s+specified.*",phenotype["trait"])
                    if m:
                        continue
                    x=next((m for m in r[rsID]["mappings"] if m["seq_region_name"]==chrom),None)
                    if not x:
                        continue
                    df.loc[i]=[rsID,r[rsID]["most_severe_consequence"].replace("_"," "),chrom+":"+str(x["start"]),phenotype["trait"],phenotype["source"],str(abs(pos - int(x["start"])))]
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

# for a list of variant mappings, return unique chr:pos pairs
def getChrPosList(mappings):
    L=list()
    for m in mappings:
        t=(m["chr"],m["pos"])
        if t not in L:
            L.append(t)

    return L

# for a list of variant mappings and a chr:pos pair, return a list of ref:alt pairs
def getMappingList(t,mappings):
    L=list()
    for m in mappings:
        if t[0]==m["chr"] and t[1]==m["pos"]:
            L.append(m)

    return L

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
    if re.search("[01]\.\d+",str(data["MAF"])):
        res["MAF"]=str(data["MAF"])
    else:
        res["MAF"]="NA"

    res["rsID"]=rs
    res["class"]=data["var_class"]
    res["consequence"]=data["most_severe_consequence"]
    if "synonyms" in data:
        res["synonyms"]=list(filter(lambda x:x!=rs,data["synonyms"]))
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

#---------------- chr:pos dependent scores -----------------

    scores=dict()
    for m in mappings:
        scores[m["chr"]+":"+str(m["pos"])]={"avg_gerp":"NA","gerp":"NA","gwava":"NA"}

#-----------------------------------------------------

    res["mappings"]=mappings
    res["population_data"]=population_data
    res["phenotype_data"]=phenotype_data
    res["clinical_significance"]=clinical_significance
    res["scores"]=scores
    
    return res

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
                        S.add(var)
                        break
    return S

# ==============================================================================================================================

# for all chr:pos mappings
def getGwavaScore(variant_data):
    chrpos=getChrPosList(variant_data["mappings"])

    for i in range(0,len(chrpos)):
        keystr=chrpos[i][0]+":"+str(chrpos[i][1])
        LOGGER.debug("GWAVA scores for mapping %s\n" %keystr)

        in_bed=tf.NamedTemporaryFile(delete=False,mode="w")
        LOGGER.debug("Input b38 bed file: %s" % in_bed.name)
        bed37_fname=tf.mktemp()
        LOGGER.debug("Output b37 bed file: %s" % bed37_fname)
        unmapped_fname=tf.mktemp()
        LOGGER.debug("Unmapped file: %s" % unmapped_fname)

        prefix=""
        if re.search("^\d+$",chrpos[i][0]) or re.search("^[XYM]$",chrpos[i][0]):
            prefix="chr"

        in_bed.write("%s%s\t%d\t%d\t%s\n" %(prefix,chrpos[i][0],chrpos[i][1]-1,chrpos[i][1],variant_data["rsID"]))
        in_bed.close()

        LOGGER.debug("Calling: liftOver %s /usr/bin/hg38ToHg19.over.chain.gz %s %s\n" %(in_bed.name,bed37_fname,unmapped_fname))
        cmdline="liftOver %s /usr/bin/hg38ToHg19.over.chain.gz %s %s" %(in_bed.name,bed37_fname,unmapped_fname)
        subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()

        if not os.path.isfile(bed37_fname):
            LOGGER.error("liftOver failed to create output file %s" % bed37_fname)
            continue

        if os.path.getsize(bed37_fname)==0:
            LOGGER.warning("liftOver produced empty output file %s" % bed37_fname)
            continue
    
        env=os.environ.copy()
        env["GWAVA_DIR"]=config.GWAVA_DIR

        annot_fname=tf.mktemp()
        LOGGER.debug("GWAVA annotation file: %s" % annot_fname)
        gwava_fname=tf.mktemp()
        LOGGER.debug("GWAVA score file: %s" % gwava_fname)

        LOGGER.debug("Calling: python2 %s/src/gwava_annotate.py %s %s" %(config.GWAVA_DIR,bed37_fname,annot_fname))
        cmdline="python2 %s/src/gwava_annotate.py %s %s" %(config.GWAVA_DIR,bed37_fname,annot_fname)
        subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env).communicate()

        LOGGER.debug("Calling: python2 %s/src/gwava.py tss %s %s" %(config.GWAVA_DIR,annot_fname,gwava_fname))
        cmdline="python2 %s/src/gwava.py tss %s %s" %(config.GWAVA_DIR,annot_fname,gwava_fname)
        subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env).communicate()

        for line in open(annot_fname, 'r'):
            if "start" in line: continue

            avg_gerp = line.split(",")[147]
            gerp = line.split(",")[148]

        for line in open(gwava_fname, 'r'):
            gwava = line.strip().split("\t")[4]

        variant_data["scores"][keystr]["avg_gerp"]=avg_gerp
        variant_data["scores"][keystr]["gerp"]=gerp
        variant_data["scores"][keystr]["gwava"]=gwava

        LOGGER.info("avg_gerp: %s: gerp: %s; gwava: %s" % (avg_gerp,gerp,gwava))

        LOGGER.debug("Removing temporary files\n")
        if os.path.isfile(in_bed.name):
            os.remove(in_bed.name)
        if os.path.isfile(bed37_fname):
            os.remove(bed37_fname)
        if os.path.isfile(unmapped_fname):
            os.remove(unmapped_fname)
        if os.path.isfile(annot_fname):
            os.remove(annot_fname)
        if os.path.isfile(gwava_fname):
            os.remove(gwava_fname)

# ==============================================================================================================================

# TODO: add links like in the original version
# TODO: gwava and gerp
def variant2df(var_data,mappings):
    df=pd.DataFrame(columns=["Value"])
    df.loc["ID"]=[var_data["rsID"]]
    keystr=mappings[0]["chr"]+":"+str(mappings[0]["pos"])
    df.loc["Location"]=[keystr]
    L=list()
    for m in mappings:
        L.append(m["ref"]+"/"+m["alt"])
    df.loc["Allele string(s)"]=[",".join(L)]
    # df.loc["SIFT score"]=[mapping["sift_score"]]
    # df.loc["SIFT prediction"]=[mapping["sift_prediction"]]
    # df.loc["PolyPhen score"]=[mapping["polyphen_score"]]
    # df.loc["PolyPhen prediction"]=[mapping["polyphen_prediction"]]
    df.loc["Average GERP score"]=var_data["scores"][keystr]["avg_gerp"]
    df.loc["GERP score"]=var_data["scores"][keystr]["gerp"]
    df.loc["GWAVA score"]=var_data["scores"][keystr]["gwava"]
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

def population2df(pop_data):
    # all alleles in all populations
    all_alleles=set()
    for z in pop_data:
        for a in z["frequency"]:
            all_alleles.add(a)

    A=list(all_alleles)
    A.sort()

    C=["Population"]
    for i in range(1,len(all_alleles)+1):
        C.append("Allele %d" %i)

    df=pd.DataFrame(columns=C)
    if len(pop_data)==0:
        return df

    i=0

    for p in config.PopulationNames:
        #LOGGER.debug("Population: %s" %p)
        x=next((z for z in pop_data if z["population"]==p),None)
        L=[p]
        if x:
            D=x["frequency"].copy()
            for a in A:
                if a not in D:
                    D[a]=0.0000
                L.append(a+" ("+str(round(D[a],4))+")")
        else:
            for j in range(1,len(all_alleles)+1):
                L.append("NA")

        df.loc[i]=L
        i+=1

    return df
