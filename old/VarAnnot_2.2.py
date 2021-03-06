#!/usr/bin/env python3

'''
importing libraries that are independent from other functions.
This first block don't have to be executed each time:
'''

import requests, sys

# libray for parsing htlm files:
from lxml import html

# library for generate python data structure from string:
import ast

# Library to execute shell commands:
import subprocess

# Importing date library:
import datetime

# LIbrary for importing options:
import argparse

# Importing the templating library that creates html files:
egg_path='/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Django-1.8.6-py2.7.egg'
sys.path.append(egg_path)

import jinja2

from variations import *
from variations_draw import *
from genes import *
from genes_draw import *

date_string = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
version = "3.0"

# Parsing command line arguments:
parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument("-i", "--id", help="Required: input variation, rsID or SNP ID + alleles: \"14:94844947_C_T\"", required=True)
input_options.add_argument("-g", "--gwava", help="Optional: perform GWAVA prediction", required=False, action='store_true')

# Extracting command line parameters:
args = parser.parse_args()
SNP_ID = args.input
GWAVA = args.GWAVA

# Let's assume we have an rsID given by the user:
filename = SNP_ID.replace(":","_")+".html"

# Submit for retrieving variant information:
print("[Info] Retrieving variant data from Ensembl ... ",file=sys.stderr)
(population_data, variant_data) = get_rsID(SNP_ID)
print("Done",file=sys.stderr)

if GWAVA:
    print("[Info] Calculating GWAVA score of the variation ...",file=sys.stderr)
    variant_data = get_GWAVA_score(variant_data)
    print("Done",file=sys.stderr)

# Retrieving variations gwas signals in the close vicinity:
print("[Info] Finding local GWAS hits around variation ... ",file=sys.stderr)
gwas_hits = get_gwas_hits_position(variant_data)
print("Done",file=sys.stderr)

# Retrieve VEP data:
print("[Info] Calculating consequences for all overlapping transcripts ... ",file=sys.stderr)
(VEP_data, variant_data) = get_VEP_data(variant_data)
print("Done",file=sys.stderr)

# Retrieving pubmed ID-s:
print("[Info] Retrieving list of publications from Pubmed ... ",file=sys.stderr)
pubmed_ids = get_pubmedIDs(variant_data["rsID"])
print("Done",file=sys.stderr)

# Return phenotypes:
print("[Info] Returning variants with phenotype annotation ... ",file=sys.stderr)
phenotype_df = get_phenotype_vars({"chr" : variant_data["chromosome"],
               "start" : variant_data["start"] - config.WINDOW,
               "end" : variant_data["start"] + config.WINDOW,
               "var" : variant_data["start"] })
print("Done",file=sys.stderr)

# Get a list of genes close to the variation:
print("[Info] Retrieving variant data from Ensembl ... ",file=sys.stderr)
gene_list = get_gene_list(variant_data)
print("Done",file=sys.stderr)

# Get population data from the Exome Aggregation Consortium:
print("[Info] Retrieving minor allele frequencies from the Exome Aggregation Consortium data ... ",file=sys.stderr)
Exac_parsed = get_ExAC_frequencies(variant_data)
print("Done",file=sys.stderr)

# Checking variation in the GTEx database:
print("[Info] Checking variation in the GTEx database ... ",file=sys.stderr)
GTEx_genes = get_GTEx_genes(variant_data)
print("Done",file=sys.stderr)

# Get population data from the Exome Aggregation Consortium:
print("[Info] Retrieving allele frequencies from UK10K data ... ",file=sys.stderr)
ukData = get_UK10K_frequencies(variant_data)
print("Done",file=sys.stderr)

# Get regulatory features:
print("[Info] Retrieving overlapping regulatory features from Ensembl... ",file=sys.stderr)
regulation = get_regulation(variant_data)
print("Done",file=sys.stderr)

# Drawing tables for the output:
print("[Info] Formatting retrieved information ... ",file=sys.stderr)
variant_table = draw_variation_table(variant_data)
regulatory_table = draw_regulation(regulation)
gwas_table = draw_gwas_position_table(gwas_hits)
VEP_table = draw_consequence_table(VEP_data)
population_table = draw_freq_table(population_data, variant_data["allele_string"])
pubmed_table = draw_pubmed_table(pubmed_ids)
uk10K_table = draw_UK10K_table(ukData)
gene_table = draw_gene_table(gene_list, SNP_ID)
exac_table = draw_ExAC_table(Exac_parsed, variant_data["reference"])
GTEx_genes_table = draw_GTEx_eQTL_table(GTEx_genes)
phenotype_table = draw_phenotpye_table_var(phenotype_df)
footer = "Generated by VarAnnoTool v.%s on: %s" %(version, date_string)
print("Done",file=sys.stderr)

# now we create the html string based on the template file and the retrieved data:
print("[Info] Saving data... ",file=sys.stderr)
variation_template = config.VAR_TEMPLATE
variation_dict = { "variation_table" : variant_table.decode('utf-8'),
                 "footer" : footer.decode('utf-8'),
                 "gwas_table": gwas_table.decode('utf-8'),
                 "consequence_table": VEP_table.decode('utf-8'),
                 "regulation" : regulatory_table.decode('utf-8'),
                 "Genomes_freq" : population_table.decode('utf-8'),
                 "ExAC_freq" : exac_table.decode('utf-8'),
                 "UK10K_freq" : uk10K_table.decode('utf-8'),
                 "genes" : gene_table.decode('utf-8'),
                 'GTExGenes' : GTEx_genes_table.decode('utf-8'),
                 "pubmed" : pubmed_table.decode('utf-8'),
                 "phenotypes" : phenotype_table.decode('utf-8')}

html = draw_html(variation_template, variation_dict)

# Saving file:
f = open(filename, 'w')
f.write(html.encode("utf8"))

print("Done",file=sys.stderr)

print("Annotating genes ... ",file=sys.stderr)
for dist in gene_list.keys():
    for gene in gene_list[dist]:
        print("\tAnnotating %s... " % gene["ID"],file=sys.stderr)
        Annotate_gene(gene["ID"])
        print("Done",file=sys.stderr)
