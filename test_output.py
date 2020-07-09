#!/usr/bin/env python3

import requests, sys
import datetime
import argparse
import jinja2
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
#print(json.dumps(VEP_data,indent=4,sort_keys=True))
#print(mapping)
LOGGER.info("Done\n")

LOGGER.info("Retrieving list of publications from PUBMED")
# TODO: also run for synonyms
# TODO: VAR_ID can be chr_pos_ref_alt
pubmed_data=getPubmed(VAR_ID)
LOGGER.info("Got %d PMID(s)\n" %(len(pubmed_data)))

# Return phenotypes:
LOGGER.info("Retreiving nearby variants with phenotype annotation")
phenotypeDF = getVariantsWithPhenotypes(mapping["chr"],mapping["pos"])
LOGGER.info("Done\n")

# Get a list of genes close to the variation:
LOGGER.info("Retrieving genes around the variant")
gene_list=getGeneList(mapping["chr"],mapping["pos"],build=build)
LOGGER.info("Got %d genes\n" %(len(gene_list)))

# Get population data from the Exome Aggregation Consortium:
LOGGER.info("Retrieving allele frequencies from ExAC")
Exac_parsed = getExacAF(mapping["chr"],mapping["pos"],mapping["ref"],mapping["alt"])
LOGGER.info("Done\n")

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

LOGGER.info("Creating variant dataframe")
variantDF = variant2df(variant_data)

LOGGER.info("Creating regulatory dataframe")
regulationDF = regulation2df(regulation)

LOGGER.info("Creating GWAS dataframe")
gwasDF = gwas2df(gwas_hits)

LOGGER.info("Creating VEP dataframe")
vepDF = vepTranscript2df(VEP_data)

LOGGER.info("Creating populations dataframe")
populationDF = population2df(variant_data["population_data"])

LOGGER.info("Creating PubMed dataframe")
#pubmedDF = pubmed2df(pubmed_data)
pubmedDF = getPubmedDF(VAR_ID,variant_data["synonyms"])

LOGGER.info("Creating gene dataframe")
geneDF = geneList2df(gene_list)

LOGGER.info("Creating ExAC dataframe")
exacDF = exac2df(Exac_parsed)

LOGGER.info("Creating GTEx dataframe")
GTEx_genesDF = gtex2df(GTEx_genes)

# ----------------------------------------------------------------------------

D=dict()

if len(variantDF):
    D["variant_table"]=variantDF.to_html(index=True,classes='utf8',table_id="common")
if len(exacDF) and len(exacDF.columns)>1:
    D["exac_table"]=exacDF.to_html(index=False,classes='utf8',table_id="common")
if len(regulationDF):
    D["regulation_table"]=regulationDF.to_html(index=False,classes='utf8',table_id="common")
if len(gwasDF):
    D["gwas_table"]=gwasDF.to_html(index=False,classes='utf8',table_id="common")
if len(vepDF):
    D["vep_table"]=vepDF.to_html(index=False,classes='utf8',table_id="common")
if len(populationDF):
    D["population_table"]=populationDF.to_html(index=False,classes='utf8',table_id="common")
if len(pubmedDF):
    D["pubmed_table"]=pubmedDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
if phenotypeDF is not None and len(phenotypeDF):
    D["phenotype_table"]=phenotypeDF.to_html(index=False,classes='utf8',table_id="common")
if len(geneDF):
    D["gene_table"]=geneDF.to_html(index=False,classes='utf8',table_id="common")
if len(GTEx_genesDF):
    D["gtex_genes_table"]=GTEx_genesDF.to_html(index=False,classes='utf8',table_id="common")

f=open("./%s.html" % VAR_ID, 'w')
f.write(generateHTML(config.VAR_TEMPLATE,D))
f.close()

# ----------------------------------------------------------------------------

# gene_ID="ENSG00000002933"

# LOGGER.info("Retreiving general information")
# info = getGeneInfo(gene_ID,build=build)

# LOGGER.info("Retreiveing cross-references")
# xrefs = getGeneXrefs(gene_ID)

# LOGGER.info("Retreiving UniProt data")
# uniprot = getUniprotData(xrefs["UniProtKB/Swiss-Prot"][0][0])

# LOGGER.info("Retreiving GWAS data")
# gwas = gene2gwas(info["name"])

# LOGGER.info("Retreiving GTEx data")
# gtex= parseGTEx(info["chromosome"],info["start"],info["end"],gene_ID)

# LOGGER.info("Retreiving mouse data\n")
# mouseDF= getMousePhenotypes(gene_ID)

# LOGGER.info("Creating general info dataframe")
# infoDF = geneInfo2df(info)

# LOGGER.info("Creating GTEx dataframe")
# gtexDF = gtex2df(gtex)
# #s=gtexDF.style.set_table_attributes('id="common"')
# #print(s.render())

# LOGGER.info("Creating GWAS dataframe")
# gwasDF = geneGwas2df(gwas)

# LOGGER.info("Creating UniProt dataframe")
# uniprotDF = uniprot2df(uniprot)

# LOGGER.info("Creating GO terms dataframe")
# goDF = goterms2df(xrefs)

# D=dict()
# if len(infoDF):
#     D["gene_table"]=infoDF.to_html(index=False,classes='utf8',table_id="common")
# if len(goDF):
#     D["go_table"]=goDF.to_html(index=False,classes='utf8',table_id="common")
# if len(uniprotDF):
#     D["uniprot_table"]=uniprotDF.to_html(index=False,classes='utf8',table_id="common")
# if len(gwasDF):
#     D["GWAS_table"]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
# if len(gtexDF):
#     D["gtexVariants_table"]=gtexDF.to_html(index=False,classes='utf8',table_id="common")
# if len(mouseDF):
#     D["mouse_table"]=mouseDF.to_html(index=False,classes='utf8',table_id="common")

# f=open("./%s.html" % gene_ID, 'w')
# f.write(generateHTML(config.GENE_TEMPLATE,D))
# f.close()

# ----------------------------------------------------------------------------

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
