import subprocess
import logging
import os
import pandas as pd
import re

from . import config

import locale
locale.setlocale(locale.LC_CTYPE,"en_US.UTF-8")

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ================================================ GWAS CATALOG ===========================================================

def gene2gwas(gene_name):
    '''
    This function performs a query for gwas signals in a local gwas datafile.
    All reported and mapped genes will give a hit.

    Output: list of dictionaries
    '''

    #LOGGER.debug("%s %s" % (config.GWAS_FILE,gene_name))
    query = "zcat %s | fgrep -iw %s" % (config.GWAS_FILE,gene_name)

    gwas_data = []
    output = subprocess.Popen(query,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    for line in output.stdout.readlines():
        fields = line.split("\t")

        gwas_data.append({
            "rsID": fields[21],
            "p-value": fields[27],
            "Consequence": fields[24],
            "Reported genes": fields[13],
            "Mapped genes": fields[14],
            "trait": fields[7],
            "PMID": fields[1],
            "URL": fields[5],
            "Author": fields[2]
        })

    if len(gwas_data)==0:
        LOGGER.info("No GWAS signas associated with %s have been found" %(gene_name))

    return gwas_data

# ----------------------------------------------------------------------------------------------------------------------

def getGwasHits(chrom,pos,window=500000):
    '''
    Retrives a list of all gwas signals within a window around a given position

    Input: chromosome, position, window size (default: 500000)
    Output: [ {"rsID","SNPID","trait","p-value","PMID","distance"}]
    '''

    start=pos-window
    if start<1:
        start=1
    end=pos+window

    L=[]
    df=pd.read_table(config.GWAS_FILE_VAR,sep="\t",header=0,compression="gzip")
    for index, row in df[(df["CHR_ID"]==chrom) & (df["CHR_POS"].astype(int)>start) & (df["CHR_POS"].astype(int)<end)].iterrows():
        rsID=row["SNPS"]
        snpid=row["CHR_ID"]+":"+str(row["CHR_POS"])
        trait=row["DISEASE/TRAIT"]
        m=re.search("^b'(.*)'$",trait)
        if m:
            trait=m.group(1)
        pval=row["P-VALUE"]
        pmid=row["PUBMEDID"]
        dist=int(row["CHR_POS"])-pos
        if abs(dist)>1000:
            dist=str(dist//1000)+"kbp"
        else:
            dist=str(dist)+"bp"
        L.append({"rsID":rsID,"SNPID":snpid,"trait":trait,"p-value":pval,"PMID":pmid,"distance":dist})

    return L

# ----------------------------------------------------------------------------------------------------------------------

def gwas2df(gwas_data):
    df=pd.DataFrame(columns=["rsID","ID","Distance","Trait","P-value","PMID"])
    i=0
    for x in gwas_data:
        df.loc[i]=[x["rsID"],x["SNPID"],x["distance"],x["trait"],x["p-value"],x["PMID"]]
        i+=1

    return df

# ----------------------------------------------------------------------------------------------------------------------

def geneGwas2df(data):
    df=pd.DataFrame(columns=["rsID","Consequence","Reported genes","Trait","URL","Author"])
    i=0
    for x in data:
        df.loc[i]=[x["rsID"],x["Consequence"],x["Reported genes"],x["trait"],"<a href='https://"+x["URL"]+"'>Link</a>",x["Author"]]
        i+=1

    return df

