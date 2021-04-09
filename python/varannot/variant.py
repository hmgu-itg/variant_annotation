import logging
import pandas as pd
import json
import re
import os
import tempfile as tf
import subprocess

from varannot import config
from varannot import query
from varannot import utils

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def rs2spdi(ID,build="38"):
    L=[]
    z=query.restQuery(query.makeRSQueryURL(ID,build=build))
    if z:
        for x in z:
            if "spdi" in x:
                spdis=x["spdi"]
                for spdi in spdis:
                    if not spdi in L:
                        L.append(spdi)

    return L

# ==============================================================================================================================

def rsList2position(L,build="38",alleles=False):
    '''
    Input: list of rsID, build (default: 38), alleles=True/False (if we need alleles as well)
    Output: a dictionary rsID --> [{"chr":c,"pos":p}, ...], or None if query fails
    '''

    D={}
    data=utils.list2string(L)
    url=query.makeRSListQueryURL(build=build)
    z=query.restQuery(url,qtype="post",data=data)
    if z:
        for x in z:
            inputID=x["input"]
            D[inputID]=[]
            spdis=x["spdi"]
            for spdi in spdis:
                h=query.parseSPDI(spdi,build=build,alleles=alleles)
                p=h["pos"]
                c=h["chr"]
                ref=h["ref"]
                alt=h["alt"]
                z=None
                if alleles:
                    z=next((x for x in D[inputID] if x["chr"]==c and x["pos"]==p and x["ref"]==ref and x["alt"]==alt),None)
                else:
                    z=next((x for x in D[inputID] if x["chr"]==c and x["pos"]==p),None)
                if not z:
                    D[inputID].append({"chr":c,"pos":p,"ref":ref,"alt":alt})
                    
        # in case some input IDs are missing in the response
        # for ID in L:
        #     if not ID in D:
        #         D[ID]=[{"chr":None,"pos":None,"ref":None,"alt":None}]
    else:
        return None

    return D

# ==============================================================================================================================

def rs2position(ID,build="38",alleles=False):
    '''
    Given rsID, return a list of dictionaries with keys "chr", "pos"
    
    Input: rsID, build (default: 38), alleles=True/False (if we need alleles as well)
    Output: a list of dictionaries with keys "chr", "pos", or None if query fails
    '''

    L=[]
    z=query.restQuery(query.makeRSQueryURL(ID,build=build))
    if z:
        for x in z:
            spdis=x["spdi"]
            for spdi in spdis:
                LOGGER.debug("SPDI: %s" % spdi)
                h=query.parseSPDI(spdi,build=build,alleles=alleles)
                p=h["pos"]
                c=h["chr"]
                ref=h["ref"]
                alt=h["alt"]
                LOGGER.debug("%s:%d:%s:%s" % (c,p,ref,alt))
                z=None
                if alleles:
                    z=next((x for x in L if x["chr"]==c and x["pos"]==p and x["ref"]==ref and x["alt"]==alt),None)
                else:
                    z=next((x for x in L if x["chr"]==c and x["pos"]==p),None)
                if not z:
                    L.append({"chr":c,"pos":p,"ref":ref,"alt":alt})
    else:
        return None

    return L

# ==============================================================================================================================

def getVariantsWithPhenotypes(chrom,pos,window=config.PHENO_WINDOW,build="38"):
    '''
    For a given genomic region, return dataframe containing variants with phenotype annotations

    Input: chromosome, position, window (default: config.PHENO_WINDOW), build (default: "38") 
    Output: pandas dataframe with columns: "ID", "Consequence", "Location", "Phenotype", "Source", "Distance"
    '''
    
    start=pos-window
    end=pos+window

    if start<1 : 
        start=1

    empty_df=pd.DataFrame(columns=["ID","Consequence","Location","Phenotype","Source","Distance"])

    if end-start>5000000:
        LOGGER.error("Maximal region size allowed: 5Mbp")
        return empty_df

    LOGGER.debug("%s:%d; window: %d" %(chrom,pos,window))
    variants=query.restQuery(query.makePhenoOverlapQueryURL(chrom,start,end,build=build),qtype="get")
    #print(json.dumps(variants,indent=4,sort_keys=True))

    if not variants:
        return empty_df

    if len(variants)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return empty_df

    rsIDs=list()
    for var in variants:
        rsIDs.append(var["id"])

    if len(rsIDs)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return empty_df
    else:
        LOGGER.info("%d variant(s) with phenotypes were found in the region %s:%d-%d" %(len(rsIDs),chrom,start,end))

    output=[]
    i=0
    df = pd.DataFrame(columns=["ID","Consequence","Location","Phenotype","Source","Link"])
    for L in utils.chunks(rsIDs,config.VARIATION_POST_MAX):
        r=query.restQuery(query.makeRSPhenotypeQueryURL(build=build),data=utils.list2string(L),qtype="post")
        if r:
            #print(json.dumps(r,indent=4,sort_keys=True))
            for rsID in r:
                for phenotype in r[rsID]["phenotypes"]:
                    m=re.search("phenotype\s+not\s+specified",phenotype["trait"])
                    if m:
                        continue

                    x=next((m for m in r[rsID]["mappings"] if m["seq_region_name"]==chrom),None)
                    if not x:
                        continue

                    link=utils.makeLink(config.ENSEMBL_PHENO_URL %rsID,"ENSEMBL")
                    if phenotype["source"] == "ClinVar":
                        link=utils.makeLink(config.CLINVAR_URL+rsID,"ClinVar")
                    elif phenotype["source"]=="NHGRI-EBI GWAS catalog":
                        link=utils.makeLink(config.NHGRI_URL+rsID,"NHGRI-EBI")
                    df.loc[i]=[rsID,r[rsID]["most_severe_consequence"].replace("_"," "),chrom+":"+str(x["start"]),phenotype["trait"],phenotype["source"],link]
                    i+=1

    return df

# ===========================================================================================================================

def getChrPosList(mappings):
    '''
    For a list of variant mappings, return a list of unique chr:pos pairs

    Input: list of mappings
    Output: list
    '''
    L=list()
    for m in mappings:
        t=(m["chr"],m["pos"])
        if t not in L:
            L.append(t)

    return L

# ===========================================================================================================================

def getMappingList(t,mappings):
    '''
    For a list of variant mappings and a chr:pos pair, return a list of mappings corresponding to chr:pos

    Input: chr:pos tuple, list of mappings
    Output: list
    '''
    L=list()
    for m in mappings:
        if t[0]==m["chr"] and t[1]==m["pos"]:
            L.append(m)

    return L

# ===========================================================================================================================

def getVariantInfo(rs,build="38"):
    '''
    For a given variant ID, return a dictionary with variant information; keys are:
    "minor_allele"
    "MAF"
    "rsID"
    "class" : variant class
    "synonyms" : list of synonym IDs
    "consequence" : most severe consequence
    "mappings" : list of mapping dictionaries with keys: "chr", "pos", "ref", "alt", "polyphen_score", "polyphen_prediction", "sift_score", "sift_preddiction"
    "population_data" : list of dictionaries "population":{"allele":"frequency"} (from phase 3 of 1KG)
    "phenotype_data" : list of dictionaries with keys "trait", "source", "risk_allele"
    "clinical_significance" : list of clinical significance terms
    "scores" : dictionary mapping "chr:pos" string to a dictionary with keys "avg_gerp", "gerp", "gwava"
    '''

    res=dict()

    # in case provided ID is not an RS
    if not utils.isRS(rs):
        t=utils.splitID(rs)
        if t:
            return {"minor_allele":None,"MAF":None,"rsID":None,"class":rs,"synonyms":[],"consequence":None,"mappings":[{"chr":t["chr"],"pos":t["pos"],"ref":t["a1"],"alt":t["a2"],"polyphen_score":"NA","polyphen_prediction":"NA","sift_score":"NA","sift_prediction":"NA"},{"chr":t["chr"],"pos":t["pos"],"ref":t["a2"],"alt":t["a1"],"polyphen_score":"NA","polyphen_prediction":"NA","sift_score":"NA","sift_prediction":"NA"}],"population_data":None,"phenotype_data":None,"clinical_significance":None,"scores":None}
        else:
            return None

#------------------- general information ---------------

    data=query.restQuery(query.makeRsPhenotypeQuery2URL(rs,build))
    #print(json.dumps(data,indent=4,sort_keys=True))

    if not data:
        return None

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

    z=query.restQuery(query.makeRSQueryURL(rs,build=build))
    if z is None:
        return None

    for x in z:
        spdis=x["spdi"]
        for spdi in spdis:
            h=query.parseSPDI(spdi,alleles=True)
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
        #scores[m["chr"]+":"+str(m["pos"])]={"avg_gerp":"NA","gerp":"NA","gwava":"NA"}
        scores[m["chr"]+":"+str(m["pos"])]={"gwava":"NA"}

#-----------------------------------------------------

    res["mappings"]=mappings
    res["population_data"]=population_data
    res["phenotype_data"]=phenotype_data
    res["clinical_significance"]=clinical_significance
    res["scores"]=scores
    
    return res

# ===========================================================================================================================

def id2rs(varid,build="38"):
    '''
    For a given variant ID (chr_pos_A1_A2), return a set of matching rs IDs

    Input: variant ID, build (default: 38)
    Output: set of rs IDs
    '''
    S=set()

    if varid.startswith("rs"):
        return varid

    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varid)
    if not m:
        LOGGER.error("%s is malformed" % varid)
        return S

    chrom=m.group(1)
    pos=int(m.group(2))
    a1=m.group(3)
    a2=m.group(4)

    batchsize=100

    if len(a1)==1 and len(a2)==1: 
        # SNP
        r=query.restQuery(query.makeOverlapVarQueryURL(chrom,pos,pos,build=build))
        if not r:
            return S

        for v in r:
            if a1 in v["alleles"] and a2 in v["alleles"]:
                S.add(v["id"])
    else:
        # in case of indels, pull all variants around the variant's position
        window=max(len(a1),len(a2))

        r=query.restQuery(query.makeOverlapVarQueryURL(chrom,pos-window,pos+window,build=build))
        if not r:
            return S

        for v in r:
            z=query.restQuery(query.makeRSQueryURL(v["id"],build=build))
            if not z:
                continue

            for x in z:
                spdis=x["spdi"]
                var=x["id"][0]
                for spdi in spdis:
                    h=query.parseSPDI(spdi,alleles=True)
                    ref=h["ref"]
                    alt=h["alt"]
                    p=h["pos"]
                    c=h["chr"]
                    LOGGER.debug("%s : %s : %s %d %s %s" % (var,spdi,c,p,ref,alt))
                    #print(spdi)
                    #print(c,p,ref,alt,sep="\t")

                    if p!=pos:
                        continue

                    if len(ref)==1 and len(alt)==1:
                        continue

                    if (ref==a1 and alt==a2) or (ref==a2 and alt==a1):
                        S.add(var)
                        break
    return S

# ===========================================================================================================================

def id2rs_mod(varid,build="38"):
    '''
    For a given variant ID (chr_pos_A1_A2), return a set of matching rs IDs

    Input: variant ID, build (default: 38)
    Output: set of rs IDs
    '''

    S=set()

    if utils.isRS(varid):
        return {varid}

    if not utils.checkID(varid):
        LOGGER.error("Variant ID %s is malformed" % varid)
        return S

    batchsize=100

    V=utils.convertVariantID(varid)
    V1=utils.convertVariantID(varid,reverse=True)
    b=utils.checkDEL(V,build=build)
    b1=utils.checkDEL(V1,build=build)
        
    window=max(len(V["del"]),len(V["ins"]))
        
    if utils.getVarType(V)=="SNP":
        r=query.restQuery(query.makeOverlapVarQueryURL(V["seq"],V["pos"],V["pos"],build=build))
        if not r:
            return S

        for v in r:
            if V["del"] in v["alleles"] and V["ins"] in v["alleles"] and v["strand"]==1 and v["start"]==v["end"]:
                S.add(v["id"])

    else:
        r=query.restQuery(query.makeOverlapVarQueryURL(V["seq"],V["pos"]-window,V["pos"]+window,build=build))
        if not r:
            return S

        LOGGER.debug("Got %d variants around %s:%d\n" %(len(r),V["seq"],V["pos"]))
        for v in r:
            LOGGER.debug("Current variant: %s" % v["id"])
            z=query.restQuery(query.makeRSQueryURL(v["id"],build=build))
            if not z:
                continue

            for x in z:
                spdis=x["spdi"]
                var=x["id"][0]
                for spdi in spdis:
                    LOGGER.debug("SPDI: %s" % spdi)
                    V2=utils.convertSPDI(spdi,build=build)
                    LOGGER.debug("V2: %s" % V2)
                    if b:
                        if utils.equivalentVariants(V,V2,build=build):
                            S.add(var)
                            break
                    if b1:
                        if utils.equivalentVariants(V1,V2,build=build):
                            S.add(var)
                            break
                    
                        
    return S

# ===========================================================================================================================

def id2rs_spdi(varid,build="38"):
    '''
    For a given variant ID (chr_pos_A1_A2), return a set of matching rs IDs
    Variant ID is converted to SPDI which is used in a variant_recoder query

    Input: variant ID, build (default: 38)
    Output: set of rs IDs
    '''

    S=set()

    if utils.isRS(varid):
        return {varid}

    if not utils.checkID(varid):
        LOGGER.error("Variant ID %s is malformed" % varid)
        return S

    # convert ID to dict 
    V=utils.convertVariantID(varid)
    V1=utils.convertVariantID(varid,reverse=True)
    b=utils.checkDEL(V,build=build)
    b1=utils.checkDEL(V1,build=build)
        
    if utils.getVarType(V)=="SNP":
        r=query.restQuery(query.makeRSQueryURL(V["seq"],V["pos"],V["pos"],build=build))
        if not r:
            return S

        for v in r:
            if V["del"] in v["alleles"] and V["ins"] in v["alleles"] and v["strand"]==1 and v["start"]==v["end"]:
                S.add(v["id"])

    else:
        r=query.restQuery(query.makeOverlapVarQueryURL(V["seq"],V["pos"]-window,V["pos"]+window,build=build))
        if not r:
            return S

        LOGGER.debug("Got %d variants around %s:%d\n" %(len(r),V["seq"],V["pos"]))
        for v in r:
            LOGGER.debug("Current variant: %s" % v["id"])
            z=query.restQuery(query.makeRSQueryURL(v["id"],build=build))
            if not z:
                continue

            for x in z:
                spdis=x["spdi"]
                var=x["id"][0]
                for spdi in spdis:
                    LOGGER.debug("SPDI: %s" % spdi)
                    V2=utils.convertSPDI(spdi,build=build)
                    LOGGER.debug("V2: %s" % V2)
                    if b:
                        if utils.equivalentVariants(V,V2,build=build):
                            S.add(var)
                            break
                    if b1:
                        if utils.equivalentVariants(V1,V2,build=build):
                            S.add(var)
                            break
                    
                        
    return S

# ==============================================================================================================================

def getGwavaScore(variant_data):
    '''
    Annotate variant data (output of getVariantInfo) with average GWAVA score

    Input: variant data
    '''

    chrpos=getChrPosList(variant_data["mappings"])

    for i in range(0,len(chrpos)):
        keystr=chrpos[i][0]+":"+str(chrpos[i][1])
        LOGGER.debug("GWAVA scores for mapping %s\n" %keystr)

        prefix=""
        if re.search("^\d+$",chrpos[i][0]) or re.search("^[XYM]$",chrpos[i][0]):
            prefix="chr"

        L=utils.runLiftOver([{"chr":prefix+chrpos[i][0],"start":chrpos[i][1]-1,"end":chrpos[i][1],"id":variant_data["rsID"]}])
        if L is None:
            continue
        
        gwava_input=tf.NamedTemporaryFile(delete=False,mode="w")
        for x in L:
            gwava_input.write("%s\t%s\t%s\t%s\n" % (x["chr"],x["start"],x["end"],x["id"]))
        gwava_input.close()

        env=os.environ.copy()
        env["GWAVA_DIR"]=config.GWAVA_DIR

        annot_fname=tf.mktemp()
        LOGGER.debug("GWAVA annotation file: %s" % annot_fname)
        gwava_output_fname=tf.mktemp()
        LOGGER.debug("GWAVA score file: %s" % gwava_output_fname)

        LOGGER.debug("Calling: python2 %s/src/gwava_annotate.py %s %s" %(config.GWAVA_DIR,gwava_input.name,annot_fname))
        cmdline="python2 %s/src/gwava_annotate.py %s %s" %(config.GWAVA_DIR,gwava_input.name,annot_fname)
        subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env).communicate()

        LOGGER.debug("Calling: python2 %s/src/gwava.py tss %s %s" %(config.GWAVA_DIR,annot_fname,gwava_output_fname))
        cmdline="python2 %s/src/gwava.py tss %s %s" %(config.GWAVA_DIR,annot_fname,gwava_output_fname)
        subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env).communicate()

        for line in open(annot_fname, 'r'):
            if "start" in line: continue
            avg_gerp=line.split(",")[147]
            gerp=line.split(",")[148]

        for line in open(gwava_output_fname, 'r'):
            gwava=line.strip().split("\t")[4]

        #variant_data["scores"][keystr]["avg_gerp"]=avg_gerp
        #variant_data["scores"][keystr]["gerp"]=gerp
        variant_data["scores"][keystr]["gwava"]=gwava

        LOGGER.info("avg_gerp: %s: gerp: %s; gwava: %s" % (avg_gerp,gerp,gwava))

        LOGGER.debug("Removing temporary files\n")
        if os.path.isfile(annot_fname):
            os.remove(annot_fname)
        if os.path.isfile(gwava_output_fname):
            os.remove(gwava_output_fname)

# ==============================================================================================================================

def variant2df(var_data,mappings):
    '''
    Transform variant data to dataframe

    Input: variant data (output of getVariantInfo), list of variant mappings
    Output: dataframe with columns "Key", "Value"
    '''

    df=pd.DataFrame(columns=["Key","Value"])

    df.loc["ID"]=["ID",var_data["rsID"]]
    
    keystr=mappings[0]["chr"]+":"+str(mappings[0]["pos"])
    df.loc["Location"]=["Location",keystr]
    
    L=list()
    for m in mappings:
        L.append(m["ref"]+"/"+m["alt"])

    df.loc["Allele string(s)"]=["Allele string(s)",",".join(L)]
    
    # df.loc["SIFT score"]=[mapping["sift_score"]]
    # df.loc["SIFT prediction"]=[mapping["sift_prediction"]]
    # df.loc["PolyPhen score"]=[mapping["polyphen_score"]]
    # df.loc["PolyPhen prediction"]=[mapping["polyphen_prediction"]]

    #df.loc["Average GERP score"]=["Average GERP score",var_data["scores"][keystr]["avg_gerp"]]
    #df.loc["GERP score"]=["GERP score",var_data["scores"][keystr]["gerp"]]
    df.loc["GWAVA score"]=["GWAVA score",var_data["scores"][keystr]["gwava"]]
    df.loc["MAF"]=["MAF",var_data["MAF"]]
    df.loc["Consequence"]=["Consequence",var_data["consequence"]]
    df.loc["Type"]=["Type",var_data["class"]]
    
    s=", ".join(var_data["synonyms"])
    if len(var_data["synonyms"])>5:
        s=",".join(var_data["synonyms"][0:5])
        s=s+"... ("+str(len(var_data["synonyms"]))+" in total)"
    df.loc["Synonyms"]=["Synonyms",s]
    
    for p in var_data["phenotype_data"]:
        df.loc[p["trait"]+" - "+p["source"]+" - "+p["risk_allele"]]=["Phenotype",p["trait"]+" - "+p["source"]+" - "+p["risk_allele"]]
        
    if len(var_data["clinical_significance"])!=0:
        df.loc["Clinical significance"]=["Clinical significance",",".join(var_data["clinical_significance"])]

    return df
    
# ==============================================================================================================================

def population2df(pop_data,ref):
    '''
    Transform variant's population data into a dataframe

    Input: population data: "population_data" dictionary of variant data (output of getVariantInfo), REF allele
    Output: dataframe
    '''

    if not pop_data:
        return None
    
    all_alleles=set()
    for z in pop_data:
        for a in z["frequency"]:
            all_alleles.add(a)

    A=list(all_alleles)
    A.sort()

    # first element should be REF
    try:
        i=A.index(ref)
        tmp=A[0]
        A[0]=ref
        A[i]=tmp
    except:
        LOGGER.error("Could not find %s in alleles" % ref)
        C=["Population"]
        for a in A:
            C.append("%s" %a)
        df=pd.DataFrame(columns=C)
        return df
        
    C=["Population"]
    for a in A:
        C.append("%s" %a)

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
                L.append(str(round(D[a],4)))
        else:
            for j in range(0,len(A)):
                L.append("NA")

        df.loc[i]=L
        i+=1

    return df

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


