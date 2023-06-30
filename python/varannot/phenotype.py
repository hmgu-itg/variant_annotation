import logging
import pandas as pd
import json

from varannot import query

LOGGER=logging.getLogger(__name__)

# ======================================================================================================================

def byGene(gene,build="38"):
    '''
    Return phenotype annotations for a given gene

    Input: gene name, build
    Output: pandas data frame with columns rs,Location,Phenotype,P-value,Link
    '''
    
    # df=pd.DataFrame(columns=["rs","Location","Phenotype","P-value","Link"])
    r=query.restQuery(query.makeGenePhenotypeQueryURL(gene,build=build,pubmed=False))
    LOGGER.debug("%s" % json.dumps(r,indent=4,sort_keys=True))
    return None

