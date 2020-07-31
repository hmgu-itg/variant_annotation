import requests
import pandas as pd
import logging
from io import StringIO

from . import config

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==============================================================================================================================

def getUniprotData(ID):
    URL = config.UNIPROT_URL
    URL += "?query=id:" + ID
    URL += "&columns=id%2Ccomment%28FUNCTION%29%2C" # Function
    URL += "comment%28SUBUNIT%29%2C" # Subunit
    URL += "comment%28DEVELOPMENTAL%20STAGE%29%2C" # Developmental stage
    URL += "comment%28TISSUE%20SPECIFICITY%29%2C" # Tissue specificity
    URL += "comment%28CATALYTIC%20ACTIVITY%29%2C" # Catalytic activity
    URL += "comment%28DISRUPTION%20PHENOTYPE%29%2C" # Disruption phenotype
    URL += "comment%28SUBCELLULAR%20LOCATION%29%2C" # Subcellular localization
    URL += "comment%28DISEASE%29%2C" # Disease
    URL += "entry%20name&format=tab" # tab delimited format returned

    r=requests.get(URL)

    data=pd.read_csv(StringIO(r.content.decode("utf-8")),sep="\t",header=0,keep_default_na=False)

    UniprotData = {
        "Name" : data["Entry name"][0],
        "Disease" : data["Involvement in disease"][0],
        "Function": data["Function [CC]"][0],
        "Entry" : data["Entry"][0],
        "Subunit" : data["Subunit structure [CC]"][0],
        "Phenotype" : data["Disruption phenotype"][0],
        "Location" : data["Subcellular location [CC]"][0],
        "Tissue" : data["Tissue specificity"][0],
        "Development" : data["Developmental stage"][0]
    }

    return UniprotData

# ======================================================================================================================

def uniprot2df(data):
    if data is None:
        return pd.DataFrame(columns=["Field","Data"])
    df=pd.DataFrame(columns=["Field","Data"])
    i=0
    for k in data:
        df.loc[i]=[k,data[k]]
        i+=1

    return df

