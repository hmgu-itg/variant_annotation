import pandas as pd
import logging
import config
import requests

from query import *

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ======================================================= PUBMED ===========================================================

# searches for publications containing given rsID and its synonyms, returns merged dataframe
def getPubmedDF(rsID, synonyms):
    '''
    Searches for publications containing given rsID and its synonyms, returns merged dataframe

    Input: rsID, list of synonyms
    Output: dataframe with columns "firstAuthor", "journal", "year", "URL", "title"

    '''

    LOGGER.debug("Creating DF for %s" % rsID)
    df=pubmed2df(getPubmed(rsID))
    if len(synonyms)==0:
        return df
    else:
        for v in synonyms:
            LOGGER.debug("Creating DF for %s" % v)
            df2=pubmed2df(getPubmed(v))
            LOGGER.debug("Merging")
            df=pd.concat([df,df2]).drop_duplicates().reset_index(drop=True)

    return df

def getPubmed(rsID):
    '''
    This function returns a list of PMIDs of those publications where the given rsID was mentioned
    Up to 1000 IDs are returned

    Input  : rsID
    Output : dictionary with PMIDs as keys and dictionaries {"firstAuthor", "title", "journal", "year", "URL"} as values
    '''
    decoded=restQuery(config.PUBMED_URL_VAR % (rsID))
    #json.dumps(decoded,indent=4,sort_keys=True)
    publication_data = {}
    if decoded is None:
        return publication_data

    pubmed_IDs=decoded["esearchresult"]["idlist"]

    for ID in pubmed_IDs:
        r=requests.get(config.PUBMED_URL_PMID % (ID))
        decoded=r.json()

        # Extracting data:
        publication_data[ID] = {
            "firstAuthor" : decoded["result"][ID]['sortfirstauthor'],
            "title" : decoded["result"][ID]['title'],
            "journal" : decoded["result"][ID]['fulljournalname'],
            "year" : decoded["result"][ID]['epubdate'].split(" ")[0],
            "URL" : "http://www.ncbi.nlm.nih.gov/pubmed/"+ID,
        }

    return publication_data

# ----------------------------------------------------------------------------------------------------------------------

def pubmed2df(pubmed_data):
    df=pd.DataFrame(columns=["First author","Journal","Year","URL","Title"])
    i=0
    for x in pubmed_data:
        d=pubmed_data[x]
        df.loc[i]=[d["firstAuthor"],d["journal"],d["year"],"<a href='"+d["URL"]+"'>Link</a>",d["title"]]
        i+=1

    return df

