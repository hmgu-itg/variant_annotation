import logging
import pandas as pd
import requests
import json

from varannot import query

LOGGER=logging.getLogger(__name__)

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
    r=requests.get(URL)
    if not r.ok:
        LOGGER.warning("Appris query failed for gene: %s (URL: %s)" % (gene_ID,URL))
        return data

    decoded=r.json()
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

def getVepDF(rsID,mappings,build="38"):
    '''
    For a given list of variant mappings (containing chr/pos/ref/alt information), 
    returns two merged dataframes
    
    Input  : rsID, list of mappings, build
    Output : dictionary ("regulatory": regulatory DataFrame,"transcript": transcript DataFrame)
    '''

    LOGGER.debug("Input: %d mapping(s)" %  len(mappings))

    output=dict()
    tr_data=pd.DataFrame(columns=["Gene ID","Transcript ID","Impact","Consequence","Isoform"])
    reg_data=pd.DataFrame(columns=["Regulatory feature ID","Impact","Consequence"])
    for m in mappings:
        data=getVepData(rsID,m,build)
        tr_data=pd.concat([tr_data,vepTranscript2df(data)]).drop_duplicates().reset_index(drop=True)
        reg_data=pd.concat([reg_data,vepRegulatory2df(data)]).drop_duplicates().reset_index(drop=True)
    output["transcript"]=tr_data
    output["regulatory"]=reg_data
    return output

# =======================================================================================================================

def getVepData(rsID,mapping_data,build="38"):
    '''
    This function returns variant's predicted effect based on chromosome, position and the alternate allele.
    Modifies input mapping_data by adding SIFT/PolyPhen scores and predictions, if available

    Input: dictionary with chr, pos, ref, alt

    Output: dictionary with keys "transcript", "regulatory"
    "transcript" : list of dictionaries with keys "ID", "impact", "gene_symbol", "gene_id", "consequence", "principal"
    "regulatory" : list of dictionaries with keys "impact", "ID", "consequence"    
    '''

    VEP_data=dict()
    VEP_data["transcript"]=[]
    VEP_data["regulatory"]=[]
    VEP=None
    
    if rsID is None:
        chrom=mapping_data["chr"]
        pos=mapping_data["pos"]
        ref=mapping_data["ref"]
        alt=mapping_data["alt"]
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
            LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt))
            return None
        VEP=query.restQuery(query.makeVepQueryURL(chrom,start,end,allele,build=build))
    else:
        VEP=query.restQuery(query.makeVepRSQueryURL(rsID,build=build))
    
    if VEP is None:
        return VEP_data

    LOGGER.debug("VEP:")
    LOGGER.debug(json.dumps(VEP,indent=4,sort_keys=True))

    sift_score=None
    polyphen_score=None
    sift_prediction=""
    polyphen_prediction=""

# ------------------------------------------------------------------------------------------------------------------

    for rec in VEP:
        if "transcript_consequences" in rec:
            for t in rec["transcript_consequences"]:
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
        
    for rec in VEP:
        if "regulatory_feature_consequences" in rec:
            for r in rec['regulatory_feature_consequences']:
                VEP_data["regulatory"].append({"impact":r["impact"],"ID":r["regulatory_feature_id"],"consequence":r["consequence_terms"]})

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

# =======================================================================================================================

def vepTranscript2df(vep_data):
    df=pd.DataFrame(columns=["Gene ID","Transcript ID","Impact","Consequence","Isoform"])
    i=0
    for x in vep_data["transcript"]:
        df.loc[i]=[x["gene_id"],x["ID"],x["impact"],x["consequence"],x["principal"]]
        i+=1

    return df

# =======================================================================================================================

def vepRegulatory2df(vep_data):
    df=pd.DataFrame(columns=["Regulatory feature ID","Impact","Consequence"])
    i=0
    for x in vep_data["regulatory"]:
        df.loc[i]=[x["ID"],x["impact"],",".join(x["consequence"])]
        i+=1

    return df

