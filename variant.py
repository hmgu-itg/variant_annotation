import config
import logging
import pandas as pd

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

#-----------------------------------------------------

    res["mappings"]=mappings
    res["population_data"]=population_data
    res["phenotype_data"]=phenotype_data
    res["clinical_significance"]=clinical_significance

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


# ==============================================================================================================================

# def get_GWAVA_score(variant_data):

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
