import pandas as pd
import logging
import json
import sys

from varannot import config
from varannot import query

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def getGeneInfo (ID,build="38"):
    '''
    Retrieve general gene information

    Input: Ensembl stable ID
    Output: dictionary
    '''
    response=query.restQuery(query.makeGeneQueryURL(ID,build=build))
    if response is None:
        return None

    gene_info=dict()
    gene_info["id"]=response["id"]
    gene_info["chromosome"]=response["seq_region_name"]
    gene_info["start"]=response["start"]
    gene_info["end"]=response["end"]
    try:
        gene_info["description"]=response["description"].split("[")[0]
    except:
        gene_info["description"]="NA"
    gene_info["name"]=response["display_name"]
    gene_info["type"]=response["biotype"]
    gene_info["strand"]=response["strand"]

    return gene_info

# ==============================================================================================================================

def getGeneXrefs (ID,build="38"):
    '''
    Retrieve cross-references

    Input  : gene ID, build(default: "38")
    Output : dictionary with keys "MIM disease","MIM gene","GO","GOSlim GOA","UniProtKB/Swiss-Prot","Human Protein Atlas","ChEMBL"
                        and lists of tuples ("primary_id","description") as values  
    '''

    LOGGER.debug(ID)

    r=query.restQuery(query.makeGeneXQueryURL(ID,build=build))
    if not r:
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

# ======================================================================================================================

def getGeneList(chrom,pos,window=config.GENE_WINDOW,build="38"):
    '''
    Return all genes overlapping a window around a position

    Input: chromosome, position, window (default: config.GENE_WINDOW), build (default: "38")
    Output: list of dictionaries with keys: 
    "start", "end", "strand", "name", "description", "biotype", "ID", "distance", "orientation" 
    '''

    LOGGER.debug("Looking for genes around %s:%d; window: %d" %(chrom,pos,window))

    end=pos+window
    start=pos-window
    if start<1: 
        start=1

    overlapping_genes=query.restQuery(query.makeGeneOverlapQueryURL(chrom,start,end,build=build))
    if overlapping_genes is None:
        return None

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

def geneList2df(gene_data):
    df=pd.DataFrame(columns=["Gene name","Gene ID","Biotype","Distance","Orientation"])
    if gene_data is None:
        return df
    i=0
    for x in gene_data:
        df.loc[i]=[x["name"],x["ID"],x["biotype"],x["distance"],x["orientation"]]
        i+=1

    return df

# ======================================================================================================================

def geneInfo2df(gene_info):
    df=pd.DataFrame(columns=["Gene name","Description","ID","Coordinates","Strand","Type"])
    if gene_info is None:
        return df
    df.loc[0]=[gene_info["name"],gene_info["description"],gene_info["id"],gene_info["chromosome"]+":"+str(gene_info["start"])+"-"+str(gene_info["end"]),gene_info["strand"],gene_info["type"],]

    return df

# ======================================================================================================================

def goterms2df(xrefs):
    df=pd.DataFrame(columns=["GO term ID","Description"])
    if xrefs is None:
        return df
    i=0
    for x in xrefs["GO"]:
        df.loc[i]=[x[0],x[1]]
        i+=1

    return df

