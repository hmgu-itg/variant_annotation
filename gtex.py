import subprocess
import pandas as pd
import logging
import config

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ===================================================== GTEx RELATED STUFF ============================================

def gtex2df(gtex_data):
    df=pd.DataFrame(columns=["Tissue","P-value","Beta (SE)","ID","Distance from TSS"])
    i=0
    for x in gtex_data:
        for z in gtex_data[x]:
            df.loc[i]=[z["tissue"],z["p-value"],z["beta"]+" ("+z["SE"]+")",x,z["dist"]]
            i+=1

    return df

def getGTExDF(mappings):
    '''
    For a given list of variant mappings (containing chr/pos/ref/alt information), 
    return a merged dataframe
    
    Input  : list of mappings
    Output : dataframe with columns: "Tissue","P-value","Beta (SE)","ID","Distance from TSS"
    '''

    LOGGER.debug("Input: %d mappings" %  len(mappings))
    df=pd.DataFrame(columns=["Tissue","P-value","Beta (SE)","ID","Distance from TSS"])

    for m in mappings:
        df=pd.concat([df,gtex2df(parseGTEx(m["chr"],m["pos"],m["pos"],m["chr"]+"_"+str(m["pos"])+"_"+m["ref"]+"_"+m["alt"]))]).drop_duplicates().reset_index(drop=True)

    return df

# given gene ID (variant ID), retreive all variant (gene) data associated with the gene (variant): tissue, p-value
def parseGTEx(chrom,start,end,ID):
    '''
    Input  : gene ID or variant ID (format: chr_pos_ref_alt)
    Output : dictionary with variant IDs (gene IDs) as keys and dictionaries with "tissue", "p-value", "beta", "SE", "dist" as values
    '''
    filename=config.GTEX_BED
    query = "tabix %s %s:%d-%d | awk -v v=%s \'BEGIN{FS=\"\t\";}$4==v{print $0;}\'" %(filename,chrom,start,end,ID)
    output = subprocess.Popen(query.strip(),universal_newlines=True,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    D=dict()
    for line in output.stdout.readlines():
        fields=line.strip().split("\t")
        data=fields[4].split(":")
        ID2=data[0]
        tissue=data[1]
        pval=str("%.4e" % float(data[2]))
        beta=str(round(float(data[3]),4))
        SE=str(round(float(data[4]),4))
        dist=data[5]
        if not ID2 in D:
            D[ID2]=[]
        D[ID2].append({"tissue" : tissue,"p-value" : pval,"beta" : beta,"SE" : SE,"dist" : dist})

    if len(D)==0:
        LOGGER.info("No eQTL signals were found for %s" %(ID))
    
    return D
