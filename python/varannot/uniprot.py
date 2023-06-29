import requests
import pandas as pd
import logging
import json
from io import StringIO

from . import config

LOGGER=logging.getLogger(__name__)
# LOGGER.setLevel(logging.DEBUG)
# ch=logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
# formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
# ch.setFormatter(formatter)
# LOGGER.addHandler(ch)

# https://www.uniprot.org/help/return_fields
# https://rest.uniprot.org/uniprotkb/P12345?fields=id,accession,cc_tissue_specificity,cc_subcellular_location
# https://string-db.org/network/9606.ENSP00000332139

# ==============================================================================================================================

UNIPROT_FIELDS=["FUNCTION","CATALYTIC ACTIVITY","SUBCELLULAR LOCATION","DISEASE","SUBUNIT","TISSUE SPECIFICITY","DISRUPTION PHENOTYPE","PATHWAY"]

def extractUniprotInfo(data):
    d=dict()
    for r in data:
        if r["commentType"]=="FUNCTION":
            d.setdefault("FUNCTION",[]).extend([x["value"] for x in r["texts"]])
        elif r["commentType"]=="CATALYTIC ACTIVITY":
            d.setdefault("CATALYTIC ACTIVITY",[]).extend([r["reaction"]["name"]])
        elif r["commentType"]=="SUBCELLULAR LOCATION":
            d.setdefault("SUBCELLULAR LOCATION",[]).extend([x["location"]["value"] for x in r["subcellularLocations"]])
        elif r["commentType"]=="DISEASE":
            d.setdefault("DISEASE",[]).extend([r["disease"]["diseaseId"]])
        elif r["commentType"]=="SUBUNIT":
            d.setdefault("SUBUNIT",[]).extend([x["value"] for x in r["texts"]])
        elif r["commentType"]=="TISSUE SPECIFICITY":
            d.setdefault("TISSUE SPECIFICITY",[]).extend([x["value"] for x in r["texts"]])
        elif r["commentType"]=="DISRUPTION PHENOTYPE":
            d.setdefault("DISRUPTION PHENOTYPE",[]).extend([x["value"] for x in r["texts"]])
        elif r["commentType"]=="PATHWAY":
            d.setdefault("PATHWAY",[]).extend([x["value"] for x in r["texts"]])
        else:
            continue
    return d
    
def getUniprotData(ID):
    URL = config.UNIPROT_URL
    # URL += "?query=id:" + ID
    URL += ID+"?fields=cc_tissue_specificity,cc_subcellular_location,cc_catalytic_activity,cc_function,cc_pathway,cc_subunit,cc_disruption_phenotype,cc_disease"
    URL += "&format=json"
    # URL += "?columns=id%2Ccomment%28FUNCTION%29%2C" # Function
    # URL += "comment%28SUBUNIT%29%2C" # Subunit
    # URL += "comment%28DEVELOPMENTAL%20STAGE%29%2C" # Developmental stage
    # URL += "comment%28TISSUE%20SPECIFICITY%29%2C" # Tissue specificity
    # URL += "comment%28CATALYTIC%20ACTIVITY%29%2C" # Catalytic activity
    # URL += "comment%28DISRUPTION%20PHENOTYPE%29%2C" # Disruption phenotype
    # URL += "comment%28SUBCELLULAR%20LOCATION%29%2C" # Subcellular localization
    # URL += "comment%28DISEASE%29%2C" # Disease
    # URL += "entry%20name&format=tsv" # tab delimited format returned

    try:
        r=requests.get(URL)
    except Exception as e:
        LOGGER.error(type(e).__name__)
        return None

    LOGGER.debug(URL)
    text=r.json()
    # LOGGER.debug(str(r.content.decode("utf-8")))    
    # data=pd.read_csv(StringIO(r.content.decode("utf-8")),sep="\t",header=0,keep_default_na=False)
    # return text
    LOGGER.debug("\n%s\n" % json.dumps(text,indent=4,sort_keys=True))
    return extractUniprotInfo(text["comments"])

# ======================================================================================================================

def uniprot2df(data):
    df=pd.DataFrame(columns=["Field","Value"])
    if data is None:
        return df
    i=0
    for k in UNIPROT_FIELDS:
        if k in data:
            df.loc[i]=[k,"; ".join(data[k])]
        else:
            df.loc[i]=[k,"NA"]
        i+=1
    return df

