import pandas as pd
import re
import logging

import query

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ========================================================== MOUSE STUFF ===========================================================

def getMouseID(human_ID,build="38"):
    '''Looking up mouse gene ID of a given human gene ID'''

    data=query.restQuery(query.makeHomologyURL(human_ID,build=build,species="mouse"))
    #print(json.dumps(data,indent=4,sort_keys=True))

    mouse_IDs = {}
    if len(data["data"])==0 or len(data["data"][0]["homologies"])==0:
        LOGGER.info("No mouse cross-reference for %s" %(human_ID))
        return mouse_IDs

    for homolog in data["data"][0]["homologies"]:
        try:
            z=query.restQuery(query.makeGeneQueryURL(homolog["target"]["id"],build=build))
            #print(json.dumps(z,indent=4,sort_keys=True))

            d=""
            m=re.search("(.*)\s+\[.*\]",z["description"])
            if m:
                d=m.group(1)
            mouse_IDs[z["id"]] = [d,z["display_name"]]
        except:
            continue

    return mouse_IDs

def getMgiID(mouse_ID,build="38"):
    '''Looking up MGI cross reference for a given mouse gene ID'''

    data=query.restQuery(query.makeGeneXQueryURL2(mouse_ID,build=build))
    #print(json.dumps(data,indent=4,sort_keys=True))

    # Now from the returned data we have to pick all possible homolgue IDs:
    for d in data:
        try:
            return d["primary_id"]
        except:
            continue

    LOGGER.info("MGI ID for %s not found" % (mouse_ID))
    return ""

def getMgiPhenotypes(MGI_ID):
    '''Returns phenotype information stored on http://www.informatics.jax.org/ '''

    URL = "http://www.informatics.jax.org/allele/report.txt?markerId=%s" % MGI_ID

    # Return all associated allele information:
    try:
        df=pd.read_csv(URL, sep="\t")

        # Drop unnecessary columns:
        df.drop([ "Allele Symbol", "Chromosome",  "Synonyms","Allele Attributes", "Transmission", "Unnamed: 10"], axis=1, inplace=True)
        df.columns = ["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease"]
        return df
    except:
        LOGGER.info("No phenotype was found for %s" %(MGI_ID))
        return pd.DataFrame(columns=["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease"])

def getMousePhenotypes(geneID,build="38"):
    '''returning mouse phenotype given human gene ID'''
    # Returning all mouse homologue IDs:
    mouse_gene_IDs=getMouseID(geneID,build=build)

    full_dataframe=pd.DataFrame(columns=["Allele ID", "Allele name", "Allele type","Phenotypes","Human disease","mouse gene ID","MGI ID","mouse gene name","mouse gene description"])

    if len(mouse_gene_IDs)==0:
        return full_dataframe

    MGI_IDs=dict()
    for mouse_gene_ID in mouse_gene_IDs:
        MGI_IDs[mouse_gene_ID] = getMgiID(mouse_gene_ID,build=build)

    # Once we have all the MGI identifiers, we retrieve all the phenotypes:
    for mouse_id, mgi_id in MGI_IDs.items():
        df=getMgiPhenotypes(mgi_id)

        # Adding extra columns for the record:
        df["mouse gene ID"] = mouse_id
        df["MGI ID"] = mgi_id
        df["mouse gene name"] = mouse_gene_IDs[mouse_id][1]
        df["mouse gene description"] = mouse_gene_IDs[mouse_id][0]
        full_dataframe = pd.merge(full_dataframe, df, how="outer")

    return full_dataframe

