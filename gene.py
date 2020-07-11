import pandas as pd
import logging

import config
from query import *

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

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

# ----------------------------------------------------------------------------------------------------------------------

def geneList2df(gene_data):
    df=pd.DataFrame(columns=["Gene name","Gene ID","Biotype","Distance","Orientation"])
    i=0
    for x in gene_data:
        df.loc[i]=[x["name"],x["ID"],x["biotype"],x["distance"],x["orientation"]]
        i+=1

    return df

# =================================== CONVERTING GENE RELATED DATA STRUCTURES TO DATAFRAMES =============================

def geneInfo2df(gene_info):
    df=pd.DataFrame(columns=["Gene name","Description","ID","Coordinates","Strand","Type"])
    df.loc[0]=[gene_info["name"],gene_info["description"],gene_info["id"],gene_info["chromosome"]+":"+str(gene_info["start"])+"-"+str(gene_info["end"]),gene_info["strand"],gene_info["type"],]

    return df

# ----------------------------------------------------------------------------------------------------------------------

def goterms2df(xrefs):
    df=pd.DataFrame(columns=["GO term ID","Description"])
    i=0
    for x in xrefs["GO"]:
        df.loc[i]=[x[0],x[1]]
        i+=1

    return df

