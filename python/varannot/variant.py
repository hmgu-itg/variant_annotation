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
    '''
    For a given rs ID, return a list of SPDI strings

    Input: rs ID
    Output: list of SPDI strings
    '''
    
    L=list()
    if not utils.isRS(ID):
        LOGGER.warning("%s is not a rs ID" %(ID))
        return L
        
    z=query.restQuery(query.makeRSQueryURL(ID,build=build))
    if z:
        LOGGER.debug("\n%s" % json.dumps(z,indent=4,sort_keys=True))
        for x in z:
            for x1 in x:
                if "spdi" in x[x1]:
                    spdis=x[x1]["spdi"]
                    for spdi in spdis:
                        if not spdi in L:
                            L.append(spdi)
    return L

# ==============================================================================================================================

def rsList2position(L,build="38",alleles=False):
    '''
    Input: list of rsIDs, build (default: 38), alleles=True/False (if we want to output REF/ALT alleles as well)
    Output: a dictionary rsID --> [{"chr":c,"pos":p,"ref":ref,"alt":alt}, ...], or None if query fails
    '''

    D=dict()
    z=query.restQuery(query.makeRSListQueryURL(build=build),qtype="post",data=utils.list2string(L))
    if z:
        LOGGER.debug("\n%s" % json.dumps(z,indent=4,sort_keys=True))
        for x in z:
            LOGGER.debug("\n%s" % json.dumps(x,indent=8,sort_keys=True))
            if x is None:
                continue
            for x1 in x:
                LOGGER.debug("Allele: %s" % x1)
                if x1 == "warnings":
                    for w in x[x1]:
                        LOGGER.warn("%s" % w)
                    continue
                inputID=x[x1]["input"]
                D[inputID]=[]
                spdis=x[x1]["spdi"]
                for spdi in spdis:
                    h=query.parseSPDI(spdi,build=build,alleles=alleles)
                    p=h["pos"]
                    c=h["chr"]
                    ref=h["ref"]
                    alt=h["alt"]
                    z1=None
                    if alleles:
                        z1=next((x for x in D[inputID] if x["chr"]==c and x["pos"]==p and x["ref"]==ref and x["alt"]==alt),None)
                    else:
                        z1=next((x for x in D[inputID] if x["chr"]==c and x["pos"]==p),None)
                    if z1 is None:
                        D[inputID].append({"chr":c,"pos":p,"ref":ref,"alt":alt})
    else:
        return None
    return D

# ==============================================================================================================================

def rs2position(ID,build="38",alleles=False):
    '''
    For a given rsID, return a list of dictionaries with keys chr,pos
    
    Input: rsID, build (default: 38), alleles=True/False (if we need alleles as well)
    Output: a list of dictionaries with keys "chr", "pos", "ref", "alt" or None if query fails
    '''

    L=[]
    z=query.restQuery(query.makeRSQueryURL(ID,build=build))
    if z:
        LOGGER.debug("\n%s" % json.dumps(z,indent=4,sort_keys=True))
        for x in z:
            for x1 in x:
                spdis=x[x1]["spdi"]
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
    Output: pandas dataframe with columns: "ID", "Consequence", "Location", "Phenotype", "Source", "Link"
    '''
    
    start=int(pos)-int(window)
    end=int(pos)+int(window)

    if start<1 : 
        start=1

    df=pd.DataFrame(columns=["ID","Consequence","Location","Phenotype","Source","Link"])

    if end-start>5000000:
        LOGGER.error("Maximal region size allowed: 5Mbp")
        return df

    LOGGER.debug("%s:%d; window: %d" %(chrom,int(pos),int(window)))
    variants=query.restQuery(query.makePhenoOverlapQueryURL(chrom,start,end,build=build),qtype="get")
    LOGGER.debug("\n%s" % json.dumps(variants,indent=4,sort_keys=True))

    if not variants:
        return df

    if len(variants)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return df

    rsIDs=list()
    for var in variants:
        rsIDs.append(var["id"])

    if len(rsIDs)==0: 
        LOGGER.info("No variants with phenotypes were found in the region %s:%d-%d" %(chrom,start,end))
        return df
    else:
        LOGGER.info("%d variant(s) with phenotypes were found in the region %s:%d-%d" %(len(rsIDs),chrom,start,end))

    i=0
    for L in utils.chunks(rsIDs,config.VARIATION_POST_MAX):
        r=query.restQuery(query.makeRSPhenotypeQueryURL(build=build),data=utils.list2string(L),qtype="post")
        if r:
            LOGGER.debug("\n%s" % json.dumps(r,indent=4,sort_keys=True))
            for rsID in r:
                for phenotype in r[rsID]["phenotypes"]:
                    m=re.search("phenotype\s+not\s+specified",phenotype["trait"])
                    if m:
                        continue

                    # default link
                    link=utils.makeLink(config.ENSEMBL_PHENO_URL % rsID,"ENSEMBL")
                    
                    if phenotype["source"] == "ClinVar":
                        link=utils.makeLink(config.CLINVAR_URL+rsID,"ClinVar")
                    elif phenotype["source"]=="NHGRI-EBI GWAS catalog":
                        link=utils.makeLink(config.NHGRI_URL+rsID,"NHGRI-EBI")
                        
                    df.loc[i]=[rsID,r[rsID]["most_severe_consequence"].replace("_"," "),r[rsID]["mappings"][0]["location"],phenotype["trait"],phenotype["source"],link]
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
    For a list of variant mappings and a chr:pos pair, return a list of mappings corresponding to a given chr:pos pair

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
    "mappings" : list of mapping dictionaries with keys: "chr", "pos", "ref", "alt", "polyphen_score", "polyphen_prediction", "sift_score", "sift_prediction"
    "population_data" : dictionary "alleles":A1/A2,"frequency":{"population name":"f1/f2"}
    "phenotype_data" : list of dictionaries with keys "trait", "source", "risk_allele"
    "clinical_significance" : list of clinical significance terms
    "scores" : dictionary mapping "chr:pos" string to a dictionary with keys "avg_gerp", "gerp", "gwava"
    '''

    res=dict()

    # in case provided ID is not an rs ID
    if not utils.isRS(rs):
        res.update({"minor_allele":None,"MAF":None,"rsID":None,"class":None,"synonyms":[],"consequence":None,"mappings":[],"population_data":None,"phenotype_data":None,"clinical_significance":None,"scores":None})
        R=utils.convertVariantID(rs)
        if utils.checkDEL(R,build=build):
            res["mappings"].append({"chr":R["seq"],"pos":R["pos"],"ref":R["del"],"alt":R["ins"],"polyphen_score":"NA","polyphen_prediction":"NA","sift_score":"NA","sift_prediction":"NA"})
        R=utils.convertVariantID(rs,reverse=True)
        if utils.checkDEL(R,build=build):
            res["mappings"].append({"chr":R["seq"],"pos":R["pos"],"ref":R["del"],"alt":R["ins"],"polyphen_score":"NA","polyphen_prediction":"NA","sift_score":"NA","sift_prediction":"NA"})
        return res

#------------------- general information ---------------

    data=query.restQuery(query.makeRsPhenotypeQuery2URL(rs,build))
    LOGGER.debug("\n%s" % json.dumps(data,indent=4,sort_keys=True))
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
    LOGGER.debug("\n%s" % json.dumps(z,indent=4,sort_keys=True))
    if z is None:
        return None
    for x in z:
        for k in x:
            spdis=x[k]["spdi"]
            for spdi in spdis:
                h=query.parseSPDI(spdi,alleles=True)
                ref=h["ref"]
                alt=h["alt"]
                p=h["pos"]
                c=h["chr"]
                mappings.append({"chr":c,"pos":p,"ref":ref,"alt":alt,"sift_score":"NA","sift_prediction":"NA","polyphen_score":"NA","polyphen_prediction":"NA"})        

#------------------ population data ----------------
#                
# only consider biallelic variants

    population_data=dict()
    all_alleles=set()
    for m in data["mappings"]:
        s=m["allele_string"]
        for a in s.split("/"):
            all_alleles.add(a)
    if len(all_alleles)==2:
        a1=all_alleles.pop()
        a2=all_alleles.pop()
        population_data["alleles"]=a1+"/"+a2
        population_data["frequency"]=dict()
        for pdata in data["populations"]:
            pop_name=pdata["population"]
            f=float(pdata["frequency"])
            if a1!=pdata["allele"]:
                f=1.0-f
            population_data["frequency"][pop_name]=str(f)+"/"+str(1.0-f)

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
        
        LOGGER.debug("\n%s" % json.dumps(r,indent=4,sort_keys=True))
        LOGGER.debug("Got %d variants around %s:%d\n" %(len(r),V["seq"],V["pos"]))
        for v in r:
            LOGGER.debug("Current variant: %s" % v["id"])
            z=query.restQuery(query.makeRSQueryURL(v["id"],build=build))
            if not z:
                continue

            LOGGER.debug("\n%s" % json.dumps(z,indent=4,sort_keys=True))
            for x in z:
                for x1 in x:
                    spdis=x[x1]["spdi"]
                    var=x[x1]["id"][0]
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

# add phenotype associations to a list of rs IDs
def addPhenotypesToRSList(rsIDs,build="38"):
    LOGGER.debug("Input rs list: %d variants" % len(rsIDs))
    R=dict()
    # exclude possible NAs first
    for L in utils.chunks(list(filter(lambda x:x!="NA",rsIDs)),config.VARIATION_POST_MAX):
        r=query.restQuery(query.makeRSPhenotypeQueryURL(build=build),data=utils.list2string(L),qtype="post")
        if not r is None:
            LOGGER.debug("\n=== phenotype query ====\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for v in r:
                if not v in rsIDs:
                    continue
                if "phenotypes" in r[v]:
                    R[v]=set([x["trait"] for x in list(filter(lambda x: not re.search("phenotype\s+not\s+specified",x["trait"]),r[v]["phenotypes"]))])
                else:
                    R[v]=set()
    for v in set(rsIDs)-(set(R.keys())-{"NA"}):
        R[v]=set()
    return R
    
# ===========================================================================================================================

# add VEP consequences to a list of variant IDs (1_12345_AC_A)
# most_severe_only==True: only output one gene:consequence pair, where gene's consequence is "most_severe_consequence"
def addConsequencesToIDList(varIDs,build="38",most_severe_only=False,gene_key="gene_id"):
    LOGGER.debug("Input ID list: %d variants" % len(varIDs))
    R=dict()
    # double check, make sure IDs have correct format
    for L in utils.chunks(list(filter(utils.checkID,varIDs)),config.VEP_POST_MAX//2):
        h={"variants":[]}
        for varid in L:
            V=utils.convertVariantID(varid)
            if utils.checkDEL(V,build=build):
                h["variants"].append(utils.variant2vep(V))
            else:
                V=utils.convertVariantID(varid,reverse=True)
                if utils.checkDEL(V,build=build):
                    h["variants"].append(utils.variant2vep(V))
        r=query.restQuery(query.makeVepListQueryURL(build=build),data=json.dumps(h),qtype="post")
        if not r is None:
            LOGGER.debug("\n======= VEP query ========\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for x in r:
                rs=x["id"]
                mcsq=x["most_severe_consequence"] if "most_severe_consequence" in x else "NA"
                H=dict()
                if "transcript_consequences" in x:
                    for g in x["transcript_consequences"]:
                        H.setdefault(g[gene_key],[]).extend(g["consequence_terms"])
                    for g in H:
                        H[g]=utils.getMostSevereConsequence(H[g])
                else:
                    H["NA"]=mcsq
                if most_severe_only is True:
                    if mcsq=="NA":
                        R[rs]={"NA":"NA"}
                    else:
                        g0="NA"
                        for g in H:
                            if H[g]==mcsq:
                                g0=g
                        R[rs]={g0:mcsq}
                else:
                    R[rs]=H
    s=set(varIDs)-(set(R.keys())-{"NA"})
    LOGGER.debug("No consequences found for %d IDs" % len(s))
    for v in s:
        R[v]={"NA":"NA"}
    return R

# ===========================================================================================================================

# add VEP consequences to a list of rs IDs
# most_severe_only==True: only output one gene:consequence pair, where gene's consequence is "most_severe_consequence"
def addConsequencesToRSList(rsIDs,build="38",most_severe_only=False,gene_key="gene_id"):
    LOGGER.debug("Input rs list: %d variants" % len(rsIDs))
    R=dict()
    # exclude possible NAs from the input list first
    for L in utils.chunks(list(filter(lambda x:x!="NA",rsIDs)),config.VEP_POST_MAX):
        r=query.restQuery(query.makeVepRSListQueryURL(build=build),data=utils.list2string(L),qtype="post")
        if not r is None:
            LOGGER.debug("\n======= VEP query ========\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for x in r:
                rs=x["id"]
                mcsq=x["most_severe_consequence"] if "most_severe_consequence" in x else "NA"
                H=dict()
                if "transcript_consequences" in x:
                    for g in x["transcript_consequences"]:
                        H.setdefault(g[gene_key],[]).extend(g["consequence_terms"])
                    for g in H:
                        H[g]=utils.getMostSevereConsequence(H[g])
                else:
                    H["NA"]=mcsq
                if most_severe_only is True:
                    if mcsq=="NA":
                        R[rs]={"NA":"NA"}
                    else:
                        g0="NA"
                        for g in H:
                            if H[g]==mcsq:
                                g0=g
                        R[rs]={g0:mcsq}
                else:
                    R[rs]=H
    s=set(rsIDs)-(set(R.keys())-{"NA"})
    LOGGER.debug("No consequences found for %d rs IDs" % len(s))
    for v in s:
        R[v]={"NA":"NA"}
    return R
            
# ===========================================================================================================================

def id2rs_list(varIDs,build="38",skip_non_rs=False,keep_all=True):
    H=dict()
    R=dict()
    # TODO: check ID validity and if it's an rsID
    # trying fast method first
    LOGGER.debug("Input variant list: %d elements" % len(varIDs))
    c=0
    t=2*len(varIDs)//config.VARIATION_POST_MAX
    if t%2:
        t=t+1
    for L in utils.chunks(varIDs,config.VARIATION_POST_MAX//2):
        L1=list()
        for x in L:
            # TODO: checks
            spdi=utils.var2spdi(utils.convertVariantID(x))
            H[spdi]=x
            L1.append(spdi)
            spdi=utils.var2spdi(utils.convertVariantID(x,reverse=True))
            H[spdi]=x
            L1.append(spdi)
        r=None
        while r is None:
            r=query.restQuery(query.makeRSListQueryURL(build=build),data=utils.list2string(L1),qtype="post")
            if r is None:
                LOGGER.debug("Retrying")
        for x1 in r:
            for x2 in x1:
                if "id" in x1[x2]:
                    v=H[x1[x2]["input"]]
                    if not v in R:
                        R[v]=set()
                        R[v].update(x1[x2]["id"])
                    else:
                        R[v].update(x1[x2]["id"])
        c+=1
        LOGGER.debug("Chunk %d (%d) done" % (c,t))
    LOGGER.debug("Found rsIDs for %d variants using fast method" % len(R.keys()))
    # slow method for unmapped
    unmapped=list(set(varIDs)-set(R.keys()))
    LOGGER.debug("Using slow method for %d variants" % len(unmapped))
    for v in unmapped:
        R[v]=id2rs_mod2(v,build)
    if skip_non_rs==True:
        LOGGER.debug("Filtering non rs IDs")
        for v in R:
            s=set(filter(utils.isRS,R[v]))
            if len(s)==0:
                R[v]={"NA"}
            else:
                R[v]=s
    if not keep_all is True:
        LOGGER.debug("Keeping only one rs ID")
        c=0
        for v in R:
            if len(R[v])>1:
                z=R[v].pop()
                R[v]={z}
                c+=1
        LOGGER.debug("Truncated %d sets" % c)
    return R

# ===========================================================================================================================

# faster version
def id2rs_mod2(varid,build="38"):
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

    V=utils.convertVariantID(varid)
    V1=utils.convertVariantID(varid,reverse=True)
        
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
        LOGGER.debug("\n%s" % json.dumps(r,indent=4,sort_keys=True))

        # only save indel IDs in L
        L=[]
        for v in r:
            if "alleles" in v and "id" in v:
                for a in v["alleles"]:
                    if a =="-" or len(a)>1:
                        L.append(v["id"])
                        break
        if len(L)==0:
            LOGGER.debug("No indels found")
            return S
        LOGGER.debug("%d indels found: %s" % (len(L),str(L)))
        
        # TODO: check if L is larger than allowed POST size
        z1=query.restQuery(query.makeRSListQueryURL(build=build),qtype="post",data=utils.list2string(L))
        LOGGER.debug("\n=======================\n%s\n==========================\n" % json.dumps(z1,indent=4,sort_keys=True))

        LOGGER.debug("---------- CHECK START ----------------\n")
        for v in z1:
            for x1 in v:
                if "spdi" in v[x1] and "id" in v[x1]:
                    var=v[x1]["id"][0]
                    spdis=v[x1]["spdi"]
                    for spdi in spdis:
                        V2=utils.convertSPDI(spdi,build=build)
                        LOGGER.debug("SPDI: %s; V2: %s" % (spdi,V2))
                        if utils.equivalentVariants(V,V2,build=build):
                            S.add(var)
                            break
                        if utils.equivalentVariants(V1,V2,build=build):
                            S.add(var)
                            break
        LOGGER.debug("----------- CHECK END -----------------\n")
    return S

# ===========================================================================================================================

def id2rs_spdi(varid,build="38"):
    '''
    THIS METHOD IS NOT RELIABLE FOR INDELS AS SEVERAL DIFFERENT SPDI NOTATIONS CAN BE EQUIVALENT

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
    # get SPDI records
    spdi=utils.var2spdi(V)
    spdi1=utils.var2spdi(V1)
    
    r=query.restQuery(query.makeRSQueryURL(spdi,build=build),quiet=True)
    # r is a list of dicts
    if not r is None:
        LOGGER.debug("Got results for %s" % (str(V)))
        LOGGER.debug("\n%s" % json.dumps(r,indent=4,sort_keys=True))
        for x1 in r:
            for x2 in x1:
                if "id" in x1[x2]:
                    for rs in x1[x2]["id"]:
                        S.add(rs)
    else:
        LOGGER.debug("No results for %s" % (str(V)))
                        
    r=query.restQuery(query.makeRSQueryURL(spdi1,build=build),quiet=True)
    if not r is None:
        LOGGER.debug("Got results for %s" % (str(V1)))
        LOGGER.debug("\n%s" % json.dumps(r,indent=4,sort_keys=True))
        for x1 in r:
            for x2 in x1:
                if "id" in x1[x2]:
                    for rs in x1[x2]["id"]:
                        S.add(rs)
    else:
        LOGGER.debug("No results for %s" % (str(V1)))

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

    Input: variant data (output of getVariantInfo), list of variant mappings (dicts of chr:pos:ref:alt:...)
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

def population2df(pop_data):
    '''
    Transform variant's population data into a dataframe

    Input: population data: dictionary "alleles":A1/A2,"frequency":{"population name":"f1/f2"}
    Output: dataframe with columns Population,A1,A2
    '''

    if not pop_data:
        return None
    df=pd.DataFrame(columns=["Population",pop_data["alleles"].split("/")[0],pop_data["alleles"].split("/")[1]])
    i=0
    for p in pop_data["frequency"]:
        df.loc[i]=[p,str(round(float(pop_data["frequency"][p].split("/")[0]),4)),str(round(float(pop_data["frequency"][p].split("/")[1]),4))]
        i+=1
    return df

