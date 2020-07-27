#!/usr/bin/env python3

import datetime
import argparse
import logging
import json
import os
import sys
import jinja2

from modules import config
from modules import utils
from modules import gnomad
from modules import mouse
from modules import variant
from modules import gene
from modules import regulation
from modules import gwas
from modules import gtex
from modules import vep
from modules import pubmed
from modules import uniprot

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

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

# ------------------------------------------------------------------------------------------------------------------------

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
logging.getLogger("uniprot").setLevel(verbosity)
logging.getLogger("mouse").setLevel(verbosity)

# ------------------------------------------------------------------------------------------------------------------------

if not utils.createDir(outdir):
    LOGGER.error("Could not create output dir %s" % outdir)
    sys.exit(1)

# ------------------------------------------------------------------------------------------------------------------------

date_string = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
version = "3.0"

# ------------------------------------------------------------------------------------------------------------------------

LOGGER.info("Retrieving variant data from Ensembl")
variant_data=variant.getVariantInfo(VAR_ID,build)

if GWAVA is not None:
    variant.getGwavaScore(variant_data)

chrpos=variant.getChrPosList(variant_data["mappings"])

LOGGER.info("Found %d chr:pos mapping(s)\n" %(len(chrpos)))
#print(json.dumps(variant_data,indent=4,sort_keys=True))

all_genes=list()
mapping_names=list()
D=dict()
for i in range(0,len(chrpos)):
    mappings=variant.getMappingList(chrpos[i],variant_data["mappings"])
    LOGGER.info("Current mapping: %s" %(chrpos[i][0]+":"+str(chrpos[i][1])))

    LOGGER.info("Creating variant dataframe")
    variantDF=variant.variant2df(variant_data,mappings)
    
    # Variants with phenotypes:
    LOGGER.info("Retreiving nearby variants with phenotype annotation")
    phenotypeDF=variant.getVariantsWithPhenotypes(chrpos[i][0],chrpos[i][1])
    
    # Get regulatory features
    LOGGER.info("Retrieving overlapping regulatory features")
    reg=regulation.getRegulation(chrpos[i][0],chrpos[i][1])
    LOGGER.info("Found %d overlapping regulatory feature(s)\n" %(len(reg)))
    LOGGER.info("Creating regulatory dataframe")
    regulationDF=regulation.regulation2df(reg)
    
    LOGGER.info("Looking for GWAS hits around the variant")
    gwas_hits=gwas.getGwasHits(chrpos[i][0],chrpos[i][1])
    LOGGER.info("Found %d GWAS hit(s)\n" %(len(gwas_hits)))
    LOGGER.info("Creating GWAS dataframe")
    gwasDF=gwas.gwas2df(gwas_hits)
    
    LOGGER.info("Creating VEP dataframe")
    vepDF=vep.getVepDF(mappings)["transcript"]
    #vepDF=vep["transcript"]
    
    LOGGER.info("Creating populations dataframe")
    populationDF=variant.population2df(variant_data["population_data"])

    LOGGER.info("Creating PubMed dataframe")
    pubmedDF=pubmed.getPubmedDF(VAR_ID,variant_data["synonyms"])

    # Get a list of genes close to the variation:
    LOGGER.info("Retrieving genes around the variant")
    gene_list=gene.getGeneList(chrpos[i][0],chrpos[i][1],build=build)
    all_genes.extend(gene_list)
    LOGGER.info("Got %d gene(s)\n" %(len(gene_list)))
    LOGGER.info("Creating gene dataframe")
    geneDF=gene.geneList2df(gene_list)
    
    LOGGER.info("Creating gnomAD dataframe")
    gnomadDF=gnomad.getPopulationAF(VAR_ID)
    
    LOGGER.info("Creating GTEx dataframe")
    GTEx_genesDF=gtex.getGTExDF(mappings)    
    LOGGER.info("Found %d eQTL(s)\n" % len(GTEx_genesDF))

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
f=open(config.OUTPUT_DIR+"/%s.html" %VAR_ID,"w")
f.write(utils.generateHTML(template_fname,D))
f.close()

# ----------------------------------------------------------------------------

D=dict()
gene_names=list()
for i in range(0,len(all_genes)):
    gene_ID=all_genes[i]["ID"]
    gene_names.append(gene_ID)

    LOGGER.info("Retreiving general information")
    info=gene.getGeneInfo(gene_ID,build=build)

    LOGGER.info("Retreiveing cross-references")
    xrefs=gene.getGeneXrefs(gene_ID)

    LOGGER.info("Retreiving UniProt data")
    if len(xrefs["UniProtKB/Swiss-Prot"])>0:
        up=uniprot.getUniprotData(xrefs["UniProtKB/Swiss-Prot"][0][0])
    else:
        up=None

    LOGGER.info("Retreiving GWAS data")
    gw=gwas.gene2gwas(info["name"])

    LOGGER.info("Retreiving GTEx data")
    gt=gtex.parseGTEx(info["chromosome"],info["start"],info["end"],gene_ID)

    LOGGER.info("Retreiving mouse data\n")
    mouseDF=mouse.getMousePhenotypes(gene_ID)

    LOGGER.info("Creating general info dataframe")
    infoDF=gene.geneInfo2df(info)

    LOGGER.info("Creating GTEx dataframe")
    gtexDF=gtex.gtex2df(gt)

    LOGGER.info("Creating GWAS dataframe")
    gwasDF=gwas.geneGwas2df(gw)

    LOGGER.info("Creating UniProt dataframe")
    uniprotDF=uniprot.uniprot2df(up)

    LOGGER.info("Creating GO terms dataframe")
    goDF=gene.goterms2df(xrefs)

    if len(infoDF)>0:
        D["gene_table%d" %i]=infoDF.to_html(index=False,classes='utf8',table_id="common")
    if len(goDF):
        D["go_table%d" %i]=goDF.to_html(index=False,classes='utf8',table_id="common")
    if uniprotDF is not None and len(uniprotDF)>0:
        D["uniprot_table%d" %i]=uniprotDF.to_html(index=False,classes='utf8',table_id="common")
    if len(gwasDF)>0:
        D["gwas_table%d" %i]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
    if len(gtexDF)>0:
        D["gtexVariants_table%d" %i]=gtexDF.to_html(index=False,classes='utf8',table_id="common")
    if len(mouseDF)>0:
        D["mouse_table%d" %i]=mouseDF.to_html(index=False,classes='utf8',table_id="common")

# ----------------------------------------------------------------------------

template_fname=config.OUTPUT_DIR+"/template_gene.html"
utils.generateGeneTemplate(gene_names,template_fname)
f = open(config.OUTPUT_DIR+"/%s_genes.html" %VAR_ID,"w")
f.write(utils.generateHTML(template_fname,D))
f.close()

