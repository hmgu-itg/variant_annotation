import sys
import json 
import re
import subprocess
import config
import pandas as pd
from io import StringIO
import os.path
import logging
import jinja2

from query import *
from utils import *
from gtex import *
from gwas import *
from mouse import *
from pubmed import *
from exac import *
from vep import *
from uniprot import *
from gene import *

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ==============================================================================================================================

# Retrieve data from the OMIM database
# def getOmimInfo (cross_refs):
#     '''
#     Required information:
#         - list of OMIM ID-s

#     Retrieved information:
#         - OMIM entry name
#         - text field name
#         - text field (that includes references as well)

#     Retrieved data structure (one dict for all queried entries):
#     dict{
#         OMIM_ID :{
#             title : OMIM_preferredTitle
#             text : {
#                 OMIM_textSectionTitle : OMIM_textSectionContent
#     }}}

#     More information how the OMIM API works see: http://omim.org/help/api
#     The API key will expire in 2016.11.11, at that time a new key have to be required.
#     '''

#     # Extracting OMIM ID from the cross refereces data:
#     MIM_IDs = [x[1] for x in cross_refs["MIM disease"]]
#     MIM_IDs.extend([x[1] for x in cross_refs["MIM gene"]])

#     # The function returns at this point if there is not MIM IDs in the crossref
#     if len(MIM_IDs) == 0:
#         return "No OMIM entry for this gene."

#     # Constructing query string:
#     URL = config.OMIM_URL

#     for ID in MIM_IDs:
#         URL += 'mimNumber=%s&' % ID

#     # Retrieving the following fields:
#     URL += 'include=text&'
#     URL += 'include=allelicVariantList&'
#     URL += 'include=referenceList&'

#     # Adding API key:
#     URL += config.OMIM_APIKEY

#     # Required format is pyton (although it is a string indeed)ormat=python
#     URL += "format=python&"

#     # Downloading page:
#     page = requests.get(URL)

#     # Reconstructingh dictionary from the returned string:
#     OMIM_entry  = ast.literal_eval(page.content)

#     # hash to fill in:
#     OMIM_data = {}

#     # Parsing hash:
#     for entry in OMIM_entry['omim']['entryList']:
#         # GEt OMIM ID:
#         ID = entry["entry"]["mimNumber"]
#         OMIM_data[ID] = {}

#         # Get OMIM name:
#         OMIM_data[ID]['title'] = entry["entry"]["titles"]['preferredTitle']

#         # Get OMIM text:
#         OMIM_data[ID]['text'] = {}
#         for fields in entry['entry']["textSectionList"]:
#             OMIM_data[ID]['text'][fields['textSection']['textSectionTitle']] = fields['textSection']['textSectionContent']

#         # now we have to parse allelic variants:
#         # print stuff['omim']['entryList'][0]['entry']['allelicVariantList'][0]['allelicVariant'].keys()
#         # ['status', 'name', 'dbSnps', 'text', 'mutations', 'number', 'alternativeNames', 'clinvarAccessions']

#         try:
#             OMIM_data[ID]['variations'] = {}
#             for variations in entry['entry']['allelicVariantList']:
#                 OMIM_data[ID]['variations'][variations['allelicVariant']['dbSnps']] = variations['allelicVariant']['text']
#         except:
#             continue

#     return OMIM_data

# ==============================================================================================================================

# Retrievig data from the EBI's gene expression atlas:
# def getGxaData(ID):
#     '''
#     This function downloads information from the EBI's gene expression atlas.
#     Input: Ensembl ID
#     Output: pandas core series with the 10 most highly expressing tissues.
#     '''
#     # the URL pointing to the gene expression atlas:
#     URL = config.GXA_URL % ID

#     # The downloaded data will be read as a pandas dataframe:
#     try:
#         df = pd.read_csv(URL, sep='\t', comment="#")
#     except:
#         return ("No gene expression information!", "No gene information.")

#     cleanExperiments = []
#     for exp in df.Experiment.tolist():
#         m = re.search('Tissues -\s+\d*(.+)', exp)
#         cleanExperiments.append(m.groups()[0].strip())

#     df.Experiment = cleanExperiments

#     IndexNames = df.Experiment.tolist()

#     # These are the main sources that we focus on:
#     sources = ["GTEx", "FANTOM5 project - adult", "FANTOM5 project - fetal"]
#     Levels = {}

#     # Extracting top 10 tissues for each sources:
#     for source in sources:
#         try:
#             fieldName = filter(lambda x : source in x, IndexNames)[0]
#             Index = IndexNames.index(fieldName)

#             # Excising the GTEx data from dataframe:
#             GTEx = df.ix[Index]

#             GTEx = GTEx.sort_values(ascending=False, na_position='last')[0:11]
#             Levels[fieldName] = []
#             for tissue, value in GTEx.iteritems():
#                 try:
#                     if float(value) > 0 : Levels[fieldName].append([tissue, value])
#                 except:
#                     pass
#         except:
#             pass

#     return (Levels, df)

# ----------------------------------------------------------------------------------------------------------------------

