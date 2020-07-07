#!/usr/bin/env python3

import requests, sys

# libray for parsing html files:
#from lxml import html

# library for generate python data structure from string:
#import ast

# Library to execute shell commands:
#import subprocess

import datetime
import argparse

# Importing the templating library that creates html files:
#egg_path='/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Django-1.8.6-py2.7.egg'
#sys.path.append(egg_path)

import jinja2

# from variations import *
# from variations_draw import *
# from genes import *
# from genes_draw import *

from functions import *
import logging

# ------------------------------------------------------------------------------------------------------------------------

# default
build="38"

# Parsing command line arguments:
parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
input_options.add_argument("--id", "-i", help="Required: input variation, rsID or SNP ID + alleles: \"14_94844947_C_T\"", required=True)
input_options.add_argument("--gwava", "-g", help="Optional: perform GWAVA prediction", required=False, action='store_true')

# Extracting command line parameters:
args = parser.parse_args()
VAR_ID = args.id
GWAVA = args.gwava

if args.build!=None:
    build=args.build

# ------------------------------------------------------------------------------------------------------------------------

LOGGER=logging.getLogger("annotator")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# ------------------------------------------------------------------------------------------------------------------------

date_string = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
version = "3.0"

filename = VAR_ID+".html"

# ------------------------------------------------------------------------------------------------------------------------

LOGGER.info("Retrieving variant data from Ensembl")
variant_data = getVariantInfo(VAR_ID,build)
mappings=set()
for m in variant_data["mappings"]:
    t=(m["chr"],m["pos"])
    if t not in mappings:
        mappings.add(t)
LOGGER.info("Found %d mapping(s)\n" %(len(mappings)))
#print(json.dumps(variant_data,indent=4,sort_keys=True))

# there can be several chr:pos mappings
# for each chr:pos mapping there can be several ref:alt pairs

# working with only one chr:pos mapping for now

mapping=variant_data["mappings"][0]

# if GWAVA:
#     print("[Info] Calculating GWAVA score of the variation ...",file=sys.stderr)
#     variant_data = get_GWAVA_score(variant_data)
#     print("Done",file=sys.stderr)

# a variant can have more than one mapping
gwas_hits=list()
LOGGER.info("Looking for GWAS hits around the variant")
gwas_hits.extend(getGwasHits(mapping["chr"],mapping["pos"]))
LOGGER.info("Found %d GWAS hit(s)\n" %(len(gwas_hits)))
#print(gwas_hits)

LOGGER.info("Retreiving consequences for all overlapping transcripts")
# for m in variant_data["mappings"]:
#     VEP_data[m["chr"]+"_"+str(m["pos"])+"_"+m["ref"]+"_"+m["alt"]]=getVepData(m)
VEP_data=getVepData(mapping)
LOGGER.info("Done\n")

LOGGER.info("Retrieving list of publications from PUBMED")
# TODO: also run for synonyms
# TODO: VAR_ID can be chr_pos_ref_alt
pubmed_data=getPubmed(VAR_ID)
LOGGER.info("Got %d PMID(s)\n" %(len(pubmed_data)))

# # Return phenotypes:
#LOGGER.info("Returning variants with phenotype annotation")
#phenotype_df = get_phenotype_vars({"chr" : variant_data["chromosome"],
#                "start" : variant_data["start"] - config.WINDOW,
#                "end" : variant_data["start"] + config.WINDOW,
#                "var" : variant_data["start"] })
# print("Done",file=sys.stderr)

# Get a list of genes close to the variation:
LOGGER.info("Retrieving genes around the variant")
gene_list=getGeneList(mapping["chr"],mapping["pos"],build=build)
LOGGER.info("Got %d genes\n" %(len(gene_list)))

# Get population data from the Exome Aggregation Consortium:
LOGGER.info("Retrieving allele frequencies from ExAC")
Exac_parsed = getExacAF(mapping["chr"],mapping["pos"],mapping["ref"],mapping["alt"])
print("Done\n",file=sys.stderr)

# Regulated genes from GTEx
LOGGER.info("Retrieving genes regulated by the variant from the GTEx dataset")
GTEx_genes=parseGTEx(mapping["chr"],mapping["pos"],mapping["pos"],mapping["chr"]+"_"+str(mapping["pos"])+"_"+mapping["ref"]+"_"+mapping["alt"])
LOGGER.info("Got %d gene(s)\n" %(len(GTEx_genes)))

# # Get population data from the UK10K
# print("[Info] Retrieving allele frequencies from UK10K data ... ",file=sys.stderr)
# ukData = get_UK10K_frequencies(variant_data)
# print("Done",file=sys.stderr)

# Get regulatory features
LOGGER.info("Retrieving overlapping regulatory features")
regulation = getRegulation(mapping["chr"],mapping["pos"])
LOGGER.info("Found %d overlapping regulatory feature(s)\n" %(len(regulation)))

#-----------------------------------------------------------------------------------------------------------------------------

# Creating dataframes

LOGGER.info("Creating variant data frame")
variantDF = variant2df(variant_data)

LOGGER.info("Creating regulatory data frame")
regulationDF = regulation2df(regulation)

LOGGER.info("Creating GWAS data frame")
gwasDF = gwas2df(gwas_hits)

LOGGER.info("Creating VEP data frame")
vepDF = vepTranscript2df(VEP_data)

LOGGER.info("Creating populations data frame")
populationDF = population2df(variant_data["population_data"])

LOGGER.info("Creating PubMed data frame")
pubmedDF = pubmed2df(pubmed_data)

# uk10K_table = draw_UK10K_table(ukData)

LOGGER.info("Creating gene data frame")
geneDF = geneList2df(gene_list)

LOGGER.info("Creating ExAC data frame")
exacDF = exac2df(Exac_parsed)

LOGGER.info("Creating GTEx data frame")
GTEx_genesDF = gtex2df(GTEx_genes)

# LOGGER.info("Creating phenotypes data frame")
# phenotype_table = draw_phenotpye_table_var(phenotype_df)
# footer = "Generated by VarAnnoTool v.%s on: %s" %(version, date_string)
# print("Done")

# # now we create the html string based on the template file and the retrieved data:
# print("[Info] Saving data... ",file=sys.stderr)
# variation_template = config.VAR_TEMPLATE
# variation_dict = { "variation_table" : variant_table.decode('utf-8'),
#                  "footer" : footer.decode('utf-8'),
#                  "gwas_table": gwas_table.decode('utf-8'),
#                  "consequence_table": VEP_table.decode('utf-8'),
#                  "regulation" : regulatory_table.decode('utf-8'),
#                  "Genomes_freq" : population_table.decode('utf-8'),
#                  "ExAC_freq" : exac_table.decode('utf-8'),
#                  "UK10K_freq" : uk10K_table.decode('utf-8'),
#                  "genes" : gene_table.decode('utf-8'),
#                  'GTExGenes' : GTEx_genes_table.decode('utf-8'),
#                  "pubmed" : pubmed_table.decode('utf-8'),
#                  "phenotypes" : phenotype_table.decode('utf-8')}

# html = draw_html(variation_template, variation_dict)

# # Saving file:
# f = open(filename, 'w')
# f.write(html.encode("utf8"))

# print("Done",file=sys.stderr)

# print("Annotating genes ... ",file=sys.stderr)
# for dist in gene_list.keys():
#     for gene in gene_list[dist]:
#         print("\tAnnotating %s... " % gene["ID"],file=sys.stderr)
#         Annotate_gene(gene["ID"])
#         print("Done",file=sys.stderr)
