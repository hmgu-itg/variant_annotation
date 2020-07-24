#!/usr/bin/env python3

import datetime
import argparse
import logging
import json
import os

import config

import utils
import gnomad
from variant import *
from gene import *
from regulation import *
from gwas import *
from gtex import *
from vep import *
from pubmed import *

# ------------------------------------------------------------------------------------------------------------------------

# default
build="38"
verbosity=logging.INFO

# Parsing command line arguments:
parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
input_options.add_argument("--id", "-i", help="Required: input variation, rsID or SNP ID + alleles: \"14_94844947_C_T\"", required=True)
input_options.add_argument("--data", "-d", help="Required: directory with GTEx, GWAS and Ensembl Regulation data", required=True)
input_options.add_argument('--output','-o', action="store",help="Optional: output directory; defaults to variant ID in the current directory")
input_options.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
input_options.add_argument("--gwava", "-g", help="Optional: path to GWAVA directory", required=False, action='store')

# Extracting command line parameters:
args = parser.parse_args()
VAR_ID = args.id
GWAVA = args.gwava
datadir=args.data

if datadir.endswith("/"):
    datadir=datadir[:-1]

config.REGULATORY_FILE=datadir+"/regulation/regulation.bed.gz"
config.GWAS_FILE_VAR=datadir+"/gwas/gwas.tsv.gz"
config.GWAS_FILE=datadir+"/gwas/gwas_full.tsv.gz"
config.GTEX_BED=datadir+"/gtex/gtex.bed.gz"

if args.build!=None:
    build=args.build

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

outdir=os.getcwd()+"/"+VAR_ID
if args.output:
    outdir=args.output
    
if outdir.endswith("/"):
    outdir=outdir[:-1]

if not utils.createDir(outdir):
    LOGGER.error("Could not create output dir %s" % outdir)
    sys.exit(1)

config.OUTPUT_DIR=outdir

if GWAVA is not None:
    if GWAVA.endswith("/"):
        GWAVA=GWAVA[:-1]
    config.GWAVA_DIR=GWAVA

# ------------------------------------------------------------------------------------------------------------------------

LOGGER=logging.getLogger("annotator")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("variant").setLevel(verbosity)
logging.getLogger("gene").setLevel(verbosity)
logging.getLogger("regulation").setLevel(verbosity)
logging.getLogger("gwas").setLevel(verbosity)
logging.getLogger("gtex").setLevel(verbosity)
logging.getLogger("pubmed").setLevel(verbosity)
logging.getLogger("vep").setLevel(verbosity)
logging.getLogger("exac").setLevel(verbosity)
logging.getLogger("gxa").setLevel(verbosity)

# ------------------------------------------------------------------------------------------------------------------------

date_string = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
version = "3.0"

filename = VAR_ID+".html"

# ------------------------------------------------------------------------------------------------------------------------

LOGGER.info("Retrieving variant data from Ensembl")
variant_data = getVariantInfo(VAR_ID,build)

if GWAVA is not None:
    getGwavaScore(variant_data)

# for t in getChrPosList(variant_data["mappings"]):
#     print("Mapping: %s:%s" %(t[0],str(t[1])))
#     for t1 in getMappingList(t,variant_data["mappings"]):
#         print("\t%s:%s" %(t1[0],t1[1]))
#     print("")

# there can be several chr:pos mappings
# for each chr:pos there can be several ref:alt pairs

chrpos=getChrPosList(variant_data["mappings"])

LOGGER.info("Found %d chr:pos mapping(s)\n" %(len(chrpos)))
#print(json.dumps(variant_data,indent=4,sort_keys=True))

mapping_names=list()
D=dict()
for i in range(0,len(chrpos)):
    mappings=getMappingList(chrpos[i],variant_data["mappings"])
    LOGGER.info("Current mapping: %s" %(chrpos[i][0]+":"+str(chrpos[i][1])))

    LOGGER.info("Creating variant dataframe")
    variantDF = variant2df(variant_data,mappings)
    LOGGER.info("Done\n")

    # Variants with phenotypes:
    LOGGER.info("Retreiving nearby variants with phenotype annotation")
    phenotypeDF = getVariantsWithPhenotypes(chrpos[i][0],chrpos[i][1])
    LOGGER.info("Done\n")

    # Get regulatory features
    LOGGER.info("Retrieving overlapping regulatory features")
    regulation = getRegulation(chrpos[i][0],chrpos[i][1])
    LOGGER.info("Found %d overlapping regulatory feature(s)\n" %(len(regulation)))
    LOGGER.info("Creating regulatory dataframe")
    regulationDF = regulation2df(regulation)
    LOGGER.info("Done\n")

    LOGGER.info("Looking for GWAS hits around the variant")
    gwas_hits=getGwasHits(chrpos[i][0],chrpos[i][1])
    LOGGER.info("Found %d GWAS hit(s)\n" %(len(gwas_hits)))
    LOGGER.info("Creating GWAS dataframe")
    gwasDF = gwas2df(gwas_hits)
    LOGGER.info("Done\n")

    LOGGER.info("Creating VEP dataframe")
    vep = getVepDF(mappings)
    vepDF=vep["transcript"]
    LOGGER.info("Done\n")

    LOGGER.info("Creating populations dataframe")
    populationDF = population2df(variant_data["population_data"])
    LOGGER.info("Done\n")

    LOGGER.info("Creating PubMed dataframe")
    pubmedDF = getPubmedDF(VAR_ID,variant_data["synonyms"])
    LOGGER.info("Done\n")

    # Get a list of genes close to the variation:
    LOGGER.info("Retrieving genes around the variant")
    gene_list=getGeneList(chrpos[i][0],chrpos[i][1],build=build)
    LOGGER.info("Got %d gene(s)\n" %(len(gene_list)))
    LOGGER.info("Creating gene dataframe")
    geneDF = geneList2df(gene_list)
    LOGGER.info("Done\n")

    LOGGER.info("Creating gnomAD dataframe")
    gnomadDF = gnomad.getPopulationAF(VAR_ID)
    LOGGER.info("Done\n")

    LOGGER.info("Creating GTEx dataframe")
    GTEx_genesDF = getGTExDF(mappings)
    LOGGER.info("Done\n")

# ----------------------------------------------------------------------------

    if len(variantDF)>0:
        D["variant_table%d" %i]=variantDF.to_html(index=True,classes='utf8',table_id="common")
    if gnomadDF is not None:
        D["gnomad_table%d" %i]=gnomadDF.to_html(index=False,classes='utf8',table_id="common")
    if len(regulationDF)>0:
        D["regulation_table%d" %i]=regulationDF.to_html(index=False,classes='utf8',table_id="common")
    if len(gwasDF)>0:
        D["gwas_table%d" %i]=gwasDF.to_html(index=False,classes='utf8',table_id="common")
    if len(vepDF)>0:
        D["vep_table%d" %i]=vepDF.to_html(index=False,classes='utf8',table_id="common")
    if len(populationDF)>0:
        D["population_table%d" %i]=populationDF.to_html(index=False,classes='utf8',table_id="common")
    if len(pubmedDF)>0:
        D["pubmed_table%d" %i]=pubmedDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
    if phenotypeDF is not None and len(phenotypeDF)>0:
        D["phenotype_table%d" %i]=phenotypeDF.to_html(index=False,classes='utf8',table_id="common")
    if len(geneDF)>0:
        D["gene_table%d" %i]=geneDF.to_html(index=False,classes='utf8',table_id="common")
    if len(GTEx_genesDF)>0:
        D["gtex_genes_table%d" %i]=GTEx_genesDF.to_html(index=False,classes='utf8',table_id="common")
    
    mapping_names.append(chrpos[i][0]+":"+str(chrpos[i][1]))

# ----------------------------------------------------------------------------

template_fname=config.OUTPUT_DIR+"/template_var.html"
utils.generateVarTemplate(mapping_names,template_fname)
f = open(config.OUTPUT_DIR+"/%s.html" %VAR_ID,"w")
f.write(generateHTML(template_fname,D))
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
