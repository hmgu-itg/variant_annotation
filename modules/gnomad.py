import logging
import pandas as pd
import re
import json

from . import config
from . import query
from . import utils

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==============================================================================================================================

def getPopulationAF(var, dataset="gnomad_r3"):
    '''
    Get allele frequencies for a given variant

    Input  : variant ID and GnomAD dataset ID (default: "gnomad_r3")
    Output : dataframe with columns "Population", "Ref", "Alt", or empty dataframe if request fails
    '''
    
    empty_df=pd.DataFrame(columns=["Population","Ref","Alt"])

    graphql=""
    if re.search("^rs\d+$",var):
        graphql="{variant (rsid: \"%s\",dataset: %s){chrom pos ref alt genome{populations{id ac an}}}}"
    elif utils.checkGnomadID(var):
        graphql="{variant (variantId: \"%s\",dataset: %s){chrom pos ref alt genome{populations{id ac an}}}}"
    else:
        LOGGER.warning("Malformed variant ID: %s" % var)
        return empty_df
        
    req_af = {"query": graphql % (var, dataset)}

    response=query.fetchGnomAD(req_af)
    #print(json.dumps(response, indent=4, sort_keys=True))    
    if response is not None:
        afs=response["data"]
        ref=afs["variant"]["ref"]
        alt=afs["variant"]["alt"]
        df=pd.DataFrame(columns=["Population",ref,alt])
        i=0
        L=list()
        for p in afs["variant"]["genome"]["populations"]:
            L.clear()
            if not re.search("MALE",p["id"],re.I):
                L.append(p["id"])
                ac=int(p["ac"])
                an=int(p["an"])
                if ac==0:
                    L.append("1.00")
                    L.append("0.00")
                elif ac==an:
                    L.append("0.00")
                    L.append("1.00")
                else:
                    L.append("%.3e" % (float(1.0)-float(p["ac"])/float(p["an"])))
                    L.append("%.3e" % (float(p["ac"])/float(p["an"])))
                    
                df.loc[i]=L
                i+=1
        return df
    else:
        LOGGER.info("No GnomAD info for %s" % var)
        return empty_df
