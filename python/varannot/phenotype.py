import logging
import pandas as pd
import json

from varannot import query

LOGGER=logging.getLogger(__name__)

# ======================================================================================================================

def byGene(gene,build="38"):
    '''
    Return phenotype annotations for a given gene
    '''
    
    # df=pd.DataFrame(columns=["rs","Location","Phenotype","P-value","Link"])
    r=query.restQuery(query.makeGenePhenotypeQueryURL(gene,build=build,pubmed=False))
    if r:
        LOGGER.debug("%s" % json.dumps(r,indent=4,sort_keys=True))
    return r

def byRegion(chrom,start,end,build="38"):
    '''
    Return phenotype annotations for a given region
    '''
    
    r=query.restQuery(query.makeOverlapPhenotypeQueryURL(chrom,start,end,build=build))
    if r:
        LOGGER.debug("%s" % json.dumps(r,indent=4,sort_keys=True))
    return r

