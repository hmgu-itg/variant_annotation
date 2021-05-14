import subprocess
import logging
import os
import pandas as pd
import re

from varannot import config
from varannot import query
from varannot import utils

# import locale
# locale.setlocale(locale.LC_CTYPE,"en_US.UTF-8")

LOGGER=logging.getLogger(__name__)

# ======================================================================================================================

def hitsByGene(gene,build="38"):
    '''
    Return GWAS catalog phenotype annotations for a given gene

    Input: gene name, build
    Output: pandas data frame with columns rs,Location,Phenotype,P-value,Link
    '''
    
    df=pd.DataFrame(columns=["rs","Location","Phenotype","P-value","Link"])
    r=query.restQuery(query.makeGenePhenotypeQueryURL(gene,build=build))
    i=0
    if r:
        for x in r:
            if x["source"]!="NHGRI-EBI GWAS catalog":
                continue
            link="NA"
            try:
                z=x["attributes"]["external_reference"]
                link=utils.makeLink("https://pubmed.ncbi.nlm.nih.gov/"+z.split(":")[1],"PubMed")
            except KeyError:
                pass
            pval="NA"
            try:
                pval=x["attributes"]["p_value"]
            except KeyError:
                pass                
            df.loc[i]=[x["Variation"],x["location"],x["description"],pval,link]
            i+=1
    return df

# ======================================================================================================================

def hitsByRegion(chrom,pos,window=config.GWAS_WINDOW,build="38"):
    '''
    Retrives a list of all GWAS catalog signals within a window around a given position

    Input: chromosome, position, window size, build
    Output: pandas data frame with columns rs,Location,Phenotype,P-value,Link
    '''
    start=int(pos)-int(window)
    end=int(pos)+int(window)
    df=pd.DataFrame(columns=["rs","Location","Phenotype","P-value","Link"])
    r=query.restQuery(query.makeOverlapPhenotypeQueryURL(chrom,start,end,build=build))
    i=0
    if r:
        for x in r:
            for p in x["phenotype_associations"]:
                if p["source"]!="NHGRI-EBI GWAS catalog":
                    continue
                link="NA"
                try:
                    z=p["attributes"]["external_reference"]
                    link=utils.makeLink("https://pubmed.ncbi.nlm.nih.gov/"+z.split(":")[1],"PubMed")
                except KeyError:
                    pass                    
                pval="NA"
                try:
                    pval=p["attributes"]["p_value"]
                except KeyError:
                    pass                
                df.loc[i]=[x["id"],p["location"],p["description"],pval,link]
                i+=1
    return df
    
# ======================================================================================================================

# # OFFLINE VERSION
# def gene2gwas(gene_name):
#     '''
#     This function performs a query for gwas signals in a local gwas datafile.
#     All reported and mapped genes will give a hit.

#     Output: list of dictionaries
#     '''

#     #LOGGER.debug("%s %s" % (config.GWAS_FILE,gene_name))
#     query = "zcat %s | fgrep -iw %s" % (config.GWAS_FILE,gene_name)

#     gwas_data = []
#     output = subprocess.Popen(query,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
#     for line in output.stdout.readlines():
#         fields = line.split("\t")

#         gwas_data.append({
#             "rsID": fields[21],
#             "p-value": fields[27],
#             "Consequence": fields[24],
#             "Reported genes": fields[13],
#             "Mapped genes": fields[14],
#             "trait": fields[7],
#             "PMID": fields[1],
#             "URL": fields[5],
#             "Author": fields[2]
#         })

#     if len(gwas_data)==0:
#         LOGGER.info("No GWAS signals associated with %s have been found" %(gene_name))

#     return gwas_data

# # ----------------------------------------------------------------------------------------------------------------------

# # OFFLINE VERSION
# def getGwasHits(chrom,pos,window=config.GWAS_WINDOW):
#     '''
#     Retrives a list of all gwas signals within a window around a given position

#     Input: chromosome, position, window size (default: 500000)
#     Output: [ {"rsID","SNPID","trait","p-value","PMID","distance"}]
#     '''

#     start=pos-window
#     if start<1:
#         start=1
#     end=pos+window

#     L=[]
#     df=pd.read_table(config.GWAS_FILE_VAR,sep="\t",header=0,compression="gzip")
#     for index, row in df[(df["CHR_ID"]==chrom) & (df["CHR_POS"].astype(int)>start) & (df["CHR_POS"].astype(int)<end)].iterrows():
#         rsID=row["SNPS"]
#         snpid=row["CHR_ID"]+":"+str(row["CHR_POS"])
#         trait=row["DISEASE/TRAIT"]
#         m=re.search("^b'(.*)'$",trait)
#         if m:
#             trait=m.group(1)
#         pval=row["P-VALUE"]
#         pmid=row["PUBMEDID"]
#         dist=int(row["CHR_POS"])-pos
#         if abs(dist)>1000:
#             dist=str(dist//1000)+"kbp"
#         else:
#             dist=str(dist)+"bp"
#         L.append({"rsID":rsID,"SNPID":snpid,"trait":trait,"p-value":pval,"PMID":pmid,"distance":dist})

#     return L

# # ----------------------------------------------------------------------------------------------------------------------

# def gwas2df(gwas_data):
#     df=pd.DataFrame(columns=["rsID","ID","Distance","Trait","P-value","PMID"])
#     i=0
#     for x in gwas_data:
#         df.loc[i]=[x["rsID"],x["SNPID"],x["distance"],x["trait"],x["p-value"],"<a href='"+config.PUBMED_URL+str(x["PMID"])+"'>"+str(x["PMID"])+"</a>"]
#         i+=1

#     return df

# # ----------------------------------------------------------------------------------------------------------------------

# def geneGwas2df(data):
#     df=pd.DataFrame(columns=["rsID","Consequence","Reported genes","Trait","URL","Author"])
#     i=0
#     for x in data:
#         df.loc[i]=[x["rsID"],x["Consequence"],x["Reported genes"],x["trait"],"<a href='https://"+x["URL"]+"'>Link</a>",x["Author"]]
#         i+=1

#     return df

