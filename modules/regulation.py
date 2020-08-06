import subprocess
import os
import logging
import pandas as pd

from . import config

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ======================================================================================================================

def getRegulation(chrom,pos,window=config.REG_WINDOW):
    '''
    This function returns all overlapping active regulatory features within window of the variation.

    Input: chromosome, position, window (default: 2000)

    Output: dictionary with all regulatory IDs as keys, dictionary("chrom","start","end","class","cells") as values

    '''

    start=pos-window
    if start<1:
        start=1
    end=pos+window

    regulatoryFile = config.REGULATORY_FILE
    FNULL = open(os.devnull, 'w')
    query="intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\n\") -b %s -sorted" % (chrom,start,end,regulatoryFile)
    output=subprocess.Popen(query.strip(),shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=FNULL)

    parsed=dict()
    for line in output.stdout.readlines():
        #print(line)
        c0,s0,e0,c1,st,en,cl,cell,act,regid=line.strip().split("\t")
        if act=="ACTIVE":
            if regid not in parsed:
                parsed[regid]={"chrom":c1,"start":int(st),"end":int(en),"class":cl,"cells":[cell]}
            else:
                if parsed[regid]["chrom"]!=c1 or parsed[regid]["start"]!=int(st) or parsed[regid]["end"]!=int(en) or parsed[regid]["class"]!=cl:
                    LOGGER.error("Regulatory ID %s has conflicting attributes: %s %d %d %s (%s %d %d %s)" % (regid,parsed[regid]["chrom"],parsed[regid]["start"],parsed[regid]["end"],parsed[regid]["class"],c1,int(st),int(en),cell))
                    FNULL.close()
                    return None
                else:
                    if cell not in parsed[regid]["cells"]:
                        parsed[regid]["cells"].append(cell)

    FNULL.close()
    return parsed

# ============================================ CONVERTING TO DATAFRAME ============================================

def regulation2df(reg_data):
    df=pd.DataFrame(columns=["ID","Class","Cell type"])
    i=0
    for r in reg_data:
        for cell in reg_data[r]["cells"]:
            df.loc[i]=[r,reg_data[r]["class"],cell]
            i+=1
    
    return df

