#!/usr/bin/env python3

import argparse
import logging
import logging.config
import json
import os
import sys
import gzip

from modules import config
from modules import utils
from modules import gnomad
from modules import mouse
from modules import variant
from modules import gene
from modules import regulation
from modules import gwas
from modules import gxa
from modules import gtex
from modules import vep
from modules import pubmed
from modules import uniprot

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

# ------------------------------------------------------------------------------------------------------------------------

# default
build="38"
verbosity=logging.INFO

# -----------------------------------
parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument('--build','-b', action="store",help="Optional: genome build; default: 38", default="38",required=False)
input_options.add_argument("--id", "-i", help="Required: input variant, rsID or variant ID: 14_94844947_C_T",required=True)
input_options.add_argument("--data", "-d", help="Required: directory with GTEx, GWAS and Ensembl Regulation data",required=True)
input_options.add_argument('--output','-o', action="store",help="Optional: output directory; default: to current directory",required=False)
input_options.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
input_options.add_argument("--gwava", "-g", help="Optional: path to GWAVA directory", required=False, action='store')
input_options.add_argument("--html", "-html", help="Optional: output HTML instead of default JSON.GZ", required=False, action='store_true')
input_options.add_argument("--reg-window", "-reg-window", help="Optional: bp window around the input variant to look for overlapping ENSEMBL regulatory elements; default: 2000", required=False, action='store',dest="reg_window")
input_options.add_argument("--pheno-window", "-pheno-window", help="Optional: bp window around the input variant to look for variants annotated with phenotypes; default: 500000", required=False, action='store',dest="pheno_window")
input_options.add_argument("--gene-window", "-gene-window", help="Optional: bp window around the input variant to look for overlapping genes; default: 1000000", required=False, action='store',dest="gene_window")
input_options.add_argument("--gwas-window", "-gwas-window", help="Optional: bp window around the input variant to look for GWAS signals; default: 500000", required=False, action='store',dest="gwas_window")

# --------------------------------------------- PARSING COMMAND LINE ------------------------------------------------------

args=parser.parse_args()
VAR_ID=args.id
GWAVA=args.gwava
out_html=args.html
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

outdir=os.getcwd()
if args.output:
    outdir=args.output
    
if outdir.endswith("/"):
    outdir=outdir[:-1]

config.OUTPUT_DIR=outdir

logfile=outdir+"/"+VAR_ID+".log"

if GWAVA is not None:
    if GWAVA.endswith("/"):
        GWAVA=GWAVA[:-1]
    config.GWAVA_DIR=GWAVA

if args.reg_window:
    config.REG_WINDOW=args.reg_window
if args.pheno_window:
    config.PHENO_WINDOW=args.pheno_window
if args.gene_window:
    config.GENE_WINDOW=args.gene_window
if args.gwas_window:
    config.GWAS_WINDOW=args.gwas_window

# ------------------------------------------------------------------------------------------------------------------------

# logging.config.dictConfig({
#     'version': 1,
#     'disable_existing_loggers': True
# })

#LOGGER=logging.getLogger("annotator")
LOGGER=logging.getLogger()
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
#ch=logging.FileHandler(logfile)
#ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

#logging.getLogger("variant").addHandler(ch)
logging.getLogger("modules.variant").setLevel(logging.DEBUG)

# logging.getLogger("gene").addHandler(ch)
# logging.getLogger("regulation").addHandler(ch)
# logging.getLogger("gwas").addHandler(ch)
# logging.getLogger("gtex").addHandler(ch)
# logging.getLogger("pubmed").addHandler(ch)
# logging.getLogger("vep").addHandler(ch)
# logging.getLogger("gxa").addHandler(ch)
# logging.getLogger("uniprot").addHandler(ch)
# logging.getLogger("utils").addHandler(ch)
# logging.getLogger("mouse").addHandler(ch)

# ------------------------------------------------------------------------------------------------------------------------

if not utils.checkFiles([config.REGULATORY_FILE,config.GWAS_FILE_VAR,config.GWAS_FILE,config.GTEX_BED]):
    LOGGER.error("Some of the data files don't exist")
    sys.exit(1)

if not utils.createDir(outdir):
    LOGGER.error("Could not create output dir %s" % outdir)
    sys.exit(1)

# ------------------------------------------------------------------------------------------------------------------------

LOGGER.info("Retrieving variant data from ENSEMBL")
variant_data=variant.getVariantInfo(VAR_ID,build)
if variant_data is None:
    LOGGER.error("Variant data could not be retreived")
    sys.exit(1)

sys.exit(0)

if GWAVA is not None:
    LOGGER.info("Getting GWAVA scores")
    variant.getGwavaScore(variant_data)

chrpos=variant.getChrPosList(variant_data["mappings"])

all_genes=list()
mapping_names=list()
D=dict()
D1=dict()
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
    temp=vep.getVepDF(mappings)
    vepDFtr=temp["transcript"]
    vepDFreg=temp["regulatory"]
    
    LOGGER.info("Creating populations dataframe")
    populationDF=variant.population2df(variant_data["population_data"],variant_data["mappings"][0]["ref"])
    populationFig=utils.df2barchart(populationDF)

    LOGGER.info("Creating PubMed dataframe")
    pubmedDF=pubmed.getPubmedDF(VAR_ID,variant_data["synonyms"])

    # Get a list of genes close to the variation:
    LOGGER.info("Retrieving genes around the variant")
    gene_list=gene.getGeneList(chrpos[i][0],chrpos[i][1],build=build)
    if gene_list:
        all_genes.extend(gene_list)
        LOGGER.info("Got %d gene(s)\n" %(len(gene_list)))
    LOGGER.info("Creating gene dataframe")
    geneDF=gene.geneList2df(gene_list)
    
    LOGGER.info("Creating gnomAD dataframe")
    gnomadDF=gnomad.getPopulationAF(VAR_ID)
    gnomadFig=utils.df2barchart(gnomadDF)
    
    LOGGER.info("Creating GTEx dataframe")
    GTEx_genesDF=gtex.getGTExDF(mappings)    
    LOGGER.info("Found %d eQTL(s)\n" % len(GTEx_genesDF))

# ----------------------------------------------------------------------------

    if out_html:
        if len(variantDF)>0:
            D["variant_table%d" %i]=variantDF.to_html(index=True,classes='utf8',table_id="common")
        # if len(gnomadDF)>0:
        #     D["gnomad_table%d" %i]=gnomadDF.to_html(index=False,classes='utf8',table_id="common")
        if len(regulationDF)>0:
            D["regulation_table%d" %i]=regulationDF.to_html(index=False,classes='utf8',table_id="common")
        if len(gwasDF)>0:
            D["gwas_table%d" %i]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
        if len(vepDFtr)>0:
            D["vep_table%d" %i]=vepDFtr.to_html(index=False,classes='utf8',table_id="common")
        if len(vepDFreg)>0:
            D["vepreg_table%d" %i]=vepDFreg.to_html(index=False,classes='utf8',table_id="common")
        # if len(populationDF)>0:
        #     D["population_table%d" %i]=populationDF.to_html(index=False,classes='utf8',table_id="common")
        if gnomadFig is not None:
            D["gnomad_table%d" %i]=gnomadFig
        if populationFig is not None:
            D["population_table%d" %i]=populationFig
        if len(pubmedDF)>0:
            D["pubmed_table%d" %i]=pubmedDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
        if len(phenotypeDF)>0:
            D["phenotype_table%d" %i]=phenotypeDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
        if len(geneDF)>0:
            D["gene_table%d" %i]=geneDF.to_html(index=False,classes='utf8',table_id="common")
        if len(GTEx_genesDF)>0:
            D["gtex_genes_table%d" %i]=GTEx_genesDF.to_html(index=False,classes='utf8',table_id="common")
    
        mapping_names.append(chrpos[i][0]+":"+str(chrpos[i][1]))
    else:
        D1["variant_table%d" %i]=variantDF.to_json(orient="records")
        D1["gnomad_table%d" %i]=gnomadDF.to_json(orient="records")
        D1["regulation_table%d" %i]=regulationDF.to_json(orient="records")
        D1["gwas_table%d" %i]=gwasDF.to_json(orient="records")
        D1["vep_table%d" %i]=vepDFtr.to_json(orient="records")
        D1["vepreg_table%d" %i]=vepDFreg.to_json(orient="records")
        D1["population_table%d" %i]=populationDF.to_json(orient="records") # <-- TODO: add barchart
        D1["pubmed_table%d" %i]=pubmedDF.to_json(orient="records")
        D1["phenotype_table%d" %i]=phenotypeDF.to_json(orient="records")
        D1["gene_table%d" %i]=geneDF.to_json(orient="records")
        D1["gtex_genes_table%d" %i]=GTEx_genesDF.to_json(orient="records")

# ----------------------------------------------------------------------------

# if out_html:
#     template_fname=config.OUTPUT_DIR+"/template_var.html"
#     utils.generateVarTemplate(mapping_names,template_fname)
#     f=open(config.OUTPUT_DIR+"/%s.html" %VAR_ID,"w")
#     f.write(utils.generateHTML(template_fname,D))
#     f.close()

# D.clear()
gene_names=list()
LOGGER.info("Total genes: %d" % len(all_genes))
for i in range(0,len(all_genes)):
    gene_ID=all_genes[i]["ID"]
    gene_names.append(gene_ID)

    LOGGER.info("Current gene: %s",gene_ID)
    LOGGER.info("Retreiving general information")
    info=gene.getGeneInfo(gene_ID,build=build)

    LOGGER.info("Retreiveing cross-references")
    xrefs=gene.getGeneXrefs(gene_ID)
    if xrefs is None:
        up=None
    else:
        LOGGER.info("Retreiving UniProt data")
        if len(xrefs["UniProtKB/Swiss-Prot"])>0:
            up=uniprot.getUniprotData(xrefs["UniProtKB/Swiss-Prot"][0][0])
        else:
            up=None

    LOGGER.info("Retreiving GWAS data")
    gw=gwas.gene2gwas(info["name"])

    LOGGER.info("Retreiving GTEx data")
    gt=gtex.parseGTEx(info["chromosome"],info["start"],info["end"],gene_ID)

    LOGGER.info("Retreiving mouse data")
    mouseDF=mouse.getMousePhenotypes(gene_ID)

    LOGGER.info("Creating general info dataframe")
    infoDF=gene.geneInfo2df(info)

    LOGGER.info("Creating GTEx dataframe")
    gtexDF=gtex.gtex2df(gt)

    LOGGER.info("Creating GWAS dataframe")
    gwasDF=gwas.geneGwas2df(gw)

    LOGGER.info("Creating GXA dataframe")
    gxaDF=gxa.getGxaDF(gene_ID)
    gxaHeatmap=gxa.df2heatmap(gxaDF)

    LOGGER.info("Creating UniProt dataframe")
    uniprotDF=uniprot.uniprot2df(up)

    LOGGER.info("Creating GO terms dataframe")
    goDF=gene.goterms2df(xrefs)
    LOGGER.info("")

# ----------------------------------------------------------------------------

    if out_html:
        if len(infoDF)>0:
            D["gene_table_%s" % gene_ID]=infoDF.to_html(index=False,classes='utf8',table_id="common")
        if len(goDF):
            D["go_table_%s" % gene_ID]=goDF.to_html(index=False,classes='utf8',table_id="common")
        if len(uniprotDF)>0:
            D["uniprot_table_%s" % gene_ID]=uniprotDF.to_html(index=False,classes='utf8',table_id="common")
        if len(gwasDF)>0:
            D["gwas_table_%s" % gene_ID]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
        D["gxa_heatmap_%s" % gene_ID]=gxaHeatmap
        if len(gtexDF)>0:
            D["gtexVariants_table_%s" % gene_ID]=gtexDF.to_html(index=False,classes='utf8',table_id="common")
        if len(mouseDF)>0:
            D["mouse_table_%s" % gene_ID]=mouseDF.to_html(index=False,classes='utf8',table_id="common")
    else:
        D1["gene2_table_%s" %gene_ID]=infoDF.to_json(orient="records")
        D1["uniprot_table_%s" %gene_ID]=uniprotDF.to_json(orient="records")
        D1["gwas2_table_%s" %gene_ID]=gwasDF.to_json(orient="records")
        D1["gxa_table_%s" %gene_ID]=gxaDF.to_json(orient="records")
        D1["gtex_variants_table_%s" %gene_ID]=gtexDF.to_json(orient="records")
        D1["mouse_table_%s" %gene_ID]=mouseDF.to_json(orient="records")
        D1["go_table_%s" %gene_ID]=goDF.to_json(orient="records")

if out_html:
    template_fname=config.OUTPUT_DIR+"/template.html"
    utils.generateTemplate(mapping_names,gene_names,template_fname)
    f = open(config.OUTPUT_DIR+"/%s.html" %VAR_ID,"w")
    f.write(utils.generateHTML(template_fname,D))
    f.close()
else:
    LOGGER.info("Saving JSON data\n")
    with gzip.GzipFile(config.OUTPUT_DIR+"/%s.json.gz" %VAR_ID,"w") as fout:
        fout.write(json.dumps(D1).encode('utf-8'))
