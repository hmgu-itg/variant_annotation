import subprocess
import pandas as pd
import logging

from varannot import config
from varannot import utils

LOGGER=logging.getLogger(__name__)

# ===============================================================================================================================

def gtex2df(gtex_data):
    df=pd.DataFrame(columns=["Tissue","P-value","Beta (SE)","ID","Distance from TSS"])
    i=0
    for x in gtex_data:
        for z in gtex_data[x]:
            df.loc[i]=[z["tissue"],z["p-value"],z["beta"]+" ("+z["SE"]+")",x,z["dist"]]
            i+=1

    LOGGER.debug("DF: %d" % len(df))
    return df

# ================================================================================================================================

# GTEx data is in b38 coordinates, if build!="38", liftOver on mappings must be preformed
def getGTExDF(mappings,build="38"):
    '''
    For a given list of variant mappings (containing chr/pos/ref/alt information), 
    return a merged dataframe
    
    Input  : list of mappings
    Output : dataframe with columns: "Tissue","P-value","Beta (SE)","ID","Distance from TSS"
    '''

    LOGGER.debug("Input: %d mappings" %  len(mappings))
    mappings38=mappings
    if build!="38":
        mappings38=utils.liftOverMappings(mappings,source_build=build)
    df=pd.DataFrame(columns=["Tissue","P-value","Beta (SE)","ID","Distance from TSS"])

    for m in mappings38:
        df=pd.concat([df,gtex2df(parseGTEx(m["chr"],m["pos"],m["pos"],m["chr"]+"_"+str(m["pos"])+"_"+m["ref"]+"_"+m["alt"]))]).drop_duplicates().reset_index(drop=True)

    return df

# ================================================================================================================================

def parseGTEx(chrom,start,end,ID):
    '''
    Input  : gene ID or variant IDs (format: chr_pos_ref_alt/rsID)
    Output : dictionary with variant IDs (gene IDs) as keys and dictionaries with "tissue", "p-value", "beta", "SE", "dist" as values
    '''
    
    #LOGGER.debug("%s %d %d %s",chrom,start,end,ID)
    query="tabix %s %s:%d-%d" %(config.GTEX_BED,chrom,start-1,end)
    output=subprocess.Popen(query.strip(),universal_newlines=True,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    D=dict()
    for line in output.stdout.readlines():
        #LOGGER.debug(line)
        fields=line.strip().split("\t")

        match=False
        x=fields[3]
        if x==ID: # gene ID
            match=True
        else:
            z=x.split("/") # variant IDs
            if len(z)==2:
                if z[0]==ID or z[1]==ID:
                    match=True
        if not match:
            continue

        data=fields[4].split(":")
        ID2=data[0]
        tissue=data[1]
        pval=str("%.4e" % float(data[2]))
        beta=str(round(float(data[3]),4))
        SE=str(round(float(data[4]),4))
        dist=data[5]

        z=ID2.split("/") # variant IDs
        if len(z)==2:
            ID2=z[0]

        if not ID2 in D:
            D[ID2]=[]
        D[ID2].append({"tissue" : tissue,"p-value" : pval,"beta" : beta,"SE" : SE,"dist" : dist})

    return D
