import pandas as pd
import logging
import subprocess
import config
import os

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==================================================== ExAC ============================================================

def getExacDF(mappings):
    '''
    For a given list of variant mappings (containing chr/pos/ref/alt information), 
    return a merged dataframe
    
    Input  : list of mappings
    Output : dataframe with allele columns
    '''

    LOGGER.debug("Input: %d mappings" %  len(mappings))
    df=pd.DataFrame(columns=["Population"])
    columns_set=False

    for m in mappings:
        x=exac2df(getExacAF(m["chr"],m["pos"],m["ref"],m["alt"]))
        df=pd.concat([df,x]).drop_duplicates().reset_index(drop=True).fillna("0 (0)")

    return df

# -----------------------------------------------------------------------------------------------------------------------

def exac2df(exac_data):
    L=["Population"]
    L1=[]
    for x in exac_data:
        L1.append(x["allele"])
    L1.sort()
    for a in L1:
        L.append(a)

    df=pd.DataFrame(columns=L)
    ExAC_pops = ['NFE','FIN','AFR','AMR','EAS','SAS','OTH','ALL']
    i=0
    for p in ExAC_pops:
        L=[p]
        for a in L1:
            x=next((z for z in exac_data if z["allele"]==a),None)
            if x:
                if p in x["populations"]:
                    c=x["populations"][p]["count"]
                    f=x["populations"][p]["frequency"]
                    L.append("%s (%s)" %(str(c),str(round(f,4))))
                else:
                    L.append("NA (NA)")
            else:
                LOGGER.error("Could not find allele %s" %(a))

        df.loc[i]=L
        i+=1

    return df

# -----------------------------------------------------------------------------------------------------------------------

# Extract Exome Aggregation Consortium (ExAC) allele frequencies
def getExacAF(chrom,pos,ref,alt):
    '''
    For a given variant, retreive allele frequencies for all ALT alleles in all populations

    Input  : chromosome, position, REF, ALT
    Output : list of dictionaries with the following keys "allele", "populations"
    "populations" value is a dictionary with "count" and "frequency" keys
    '''

    pops=["AFR", "ALL", "FIN", "AMR", "EAS", "SAS", "OTH", "NFE"]

    ExAC_file=config.EXAC_FILE
    if not os.path.isfile(ExAC_file):
        LOGGER.error("ExAC file %s not found" %(ExAC_file))
        return None

    query="tabix %s %s:%s-%s" %(ExAC_file,chrom,str(pos),str(pos))
    output=subprocess.Popen(query.strip(),shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    parsed=list()
    for line in output.stdout.readlines():
        fields=line.strip().split("\t")
        alts = fields[4].split(",")
        if (fields[0]==chrom and int(fields[1])==pos and ref in fields[3].split(",") and alt in alts):
            data=dict()
            for x in fields[7].split(";"):
                try:
                    (key, value)=x.split("=")
                    data[key]=value.split(",")
                except:
                    pass

            for i in range(0,len(alts)):
                d=dict()
                d["allele"]=alts[i]
                d["populations"]=dict()
                for p in pops:
                    label1="AC_"+p
                    label2="AN_"+p
                    if p=="ALL": 
                        label1="AC"
                        label2="AN"
                    count="NA"
                    freq="NA"
                    if label1 in data and label2 in data:
                        count=int(data[label1][i])
                        freq=float(count)/int(data[label2][0])
                    d["populations"][p]={"count":count,"frequency":freq}
                parsed.append(d)
    
    return parsed

