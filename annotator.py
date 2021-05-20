#!/usr/bin/env python3

import argparse
import logging
import json
import os
import sys
import gzip
import pandas as pd

from varannot import config
from varannot import utils
from varannot import gnomad
from varannot import mouse
from varannot import variant
from varannot import gene
from varannot import regulation
from varannot import gwas
from varannot import gxa
from varannot import gtex
from varannot import vep
from varannot import pubmed
from varannot import uniprot

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

# ------------------------------------------------------------------------------------------------------------------------

# default
build="38"
verbosity=logging.INFO

# ------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
input_options = parser.add_argument_group('Input options')
input_options.add_argument('--build','-b', action="store",type=str,help="Optional: genome build; default: 38", default=build,required=False)
input_options.add_argument("--id", "-i", help="Required: input variant, rsID or variant ID: 14_94844947_C_T",required=True)
input_options.add_argument("--data", "-d", help="Required: directory with GWAVA, GTEx, GWAS and Ensembl Regulation data",required=True)
input_options.add_argument('--output','-o', action="store",help="Optional: output directory; default: to current directory",required=False)
input_options.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
input_options.add_argument("--json", "-j", help="Optional: output JSON.GZ instead of default HTML", required=False, action='store_true')
input_options.add_argument("--reg-window", "-reg-window", help="Optional: bp window around the input variant to look for overlapping ENSEMBL regulatory elements; default: %d" % config.REG_WINDOW, type=int, default=config.REG_WINDOW, required=False, action='store',dest="reg_window")
input_options.add_argument("--pheno-window", "-pheno-window", help="Optional: bp window around the input variant to look for variants annotated with phenotypes; default: %d" % config.PHENO_WINDOW, type=int, default=config.PHENO_WINDOW, required=False, action='store',dest="pheno_window")
input_options.add_argument("--gene-window", "-gene-window", help="Optional: bp window around the input variant to look for overlapping genes; default: %d" % config.GENE_WINDOW, type=int, default=config.GENE_WINDOW, required=False, action='store',dest="gene_window")
input_options.add_argument("--gwas-window", "-gwas-window", help="Optional: bp window around the input variant to look for GWAS signals; default: %d" % config.GWAS_WINDOW, type=int, default=config.GWAS_WINDOW, required=False, action='store',dest="gwas_window")

# --------------------------------------------- PARSING COMMAND LINE ------------------------------------------------------

args=parser.parse_args()
VAR_ID=args.id
out_html=not args.json
datadir=args.data

if not datadir.endswith("/"):
    datadir+="/"

config.REGULATORY_FILE=datadir+config.REGULATORY_FILE_SFX
config.GWAS_FILE_VAR=datadir+config.GWAS_FILE_VAR_SFX
config.GWAS_FILE=datadir+config.GWAS_FILE_SFX
config.GTEX_BED=datadir+config.GTEX_BED_SFX
config.GWAVA_DIR=datadir+"gwava"
config.GXA_FILE=datadir+config.GXA_TSV_SFX

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
config.OUTPUT_DIR_FIG=config.OUTPUT_DIR+"/figures"

logfile=outdir+"/"+VAR_ID+".log"

config.REG_WINDOW=args.reg_window
config.PHENO_WINDOW=args.pheno_window
config.GENE_WINDOW=args.gene_window
config.GWAS_WINDOW=args.gwas_window

# ------------------------------------------------------------------------------------------------------------------------

if not utils.createDir(outdir):
    print("ERROR: Could not create output dir %s" % outdir,file=sys.stderr)
    sys.exit(1)

if not utils.createDir(config.OUTPUT_DIR_FIG):
    print("ERROR: Could not create %s" % config.OUTPUT_DIR_FIG,file=sys.stderr)
    sys.exit(1)

# -------------------------------------------------------- LOGGING -------------------------------------------------------

LOGGER=logging.getLogger("annotator")
LOGGER.setLevel(verbosity)
#ch=logging.StreamHandler()
ch=logging.FileHandler(logfile,'w')
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)
LOGGER.addHandler(logging.StreamHandler(sys.stdout))

modules=["variant","gene","regulation","gwas","gtex","pubmed","uniprot","utils","mouse","vep","gnomad","query","gxa"]
for m in modules:
    logging.getLogger("varannot."+m).addHandler(ch)
    logging.getLogger("varannot."+m).addHandler(logging.StreamHandler(sys.stdout))
    logging.getLogger("varannot."+m).setLevel(verbosity)

# ------------------------------------------------------------------------------------------------------------------------

for arg in vars(args):
    LOGGER.info("INPUT OPTIONS: %s : %s" % (arg, getattr(args, arg)))
LOGGER.info("")

sys.stderr=open(logfile,"a")

# ------------------------------------------------------------------------------------------------------------------------

if not utils.checkID(VAR_ID):
    LOGGER.error("%s is not a valid variant ID" %(VAR_ID))
    sys.exit(1)

if not utils.isRS(VAR_ID) and not utils.checkAlleles(VAR_ID,build):
    LOGGER.error("Alleles in %s don't match reference sequence" %(VAR_ID))
    sys.exit(1)

if not utils.checkFiles([config.REGULATORY_FILE,config.GWAS_FILE_VAR,config.GWAS_FILE,config.GTEX_BED,config.GXA_FILE]):
    sys.exit(1)

LOGGER.info("Reading in GXA file")
config.GXA_DF=pd.read_table(config.GXA_FILE,header=0,compression="gzip")
    
# --------------------------------------------------------- MAIN ---------------------------------------------------------

LOGGER.info("Retreiving rsIDs for %s" % VAR_ID)
rsIDs=variant.id2rs_mod(VAR_ID,build=build)
LOGGER.info("Got %d rsID(s)" % len(rsIDs))
if len(rsIDs)!=0:
    LOGGER.info("rs ID(s): %s",", ".join(rsIDs))

if len(rsIDs)==0:
    rsIDs={VAR_ID}

for current_variant_id in rsIDs:
    LOGGER.info("Retrieving variant data for %s from ENSEMBL" % current_variant_id)
    variant_data=variant.getVariantInfo(current_variant_id,build)
    if variant_data is None:
        LOGGER.error("Variant data for %s could not be retreived" % current_variant_id)
        continue

    LOGGER.info("Getting GWAVA scores")
    variant.getGwavaScore(variant_data)

    # variant can have multiple chr:pos mappings
    # each chr:pos mapping can have multiple ref:alt pairs

    # all chr:pos pairs
    chrpos=variant.getChrPosList(variant_data["mappings"])

    all_genes=list()
    mapping_names=list()
    D=dict()
    D1=dict()

#----------------------------------------------------------------------------------------------------------------------

    # loop over all chr:pos mappings
    for i in range(0,len(chrpos)):
        LOGGER.info("Current chr:pos: %s" %(chrpos[i][0]+":"+str(chrpos[i][1])))
        # all chr:pos:ref:alt etc. information corresponding to the current chr:pos
        mappings=variant.getMappingList(chrpos[i],variant_data["mappings"])
        
        LOGGER.info("Creating variant dataframe")
        variantDF=variant.variant2df(variant_data,mappings)
        
        LOGGER.info("Nearby variants with phenotype annotation")
        phenotypeDF=variant.getVariantsWithPhenotypes(chrpos[i][0],chrpos[i][1])
        
        LOGGER.info("Overlapping regulatory features")
        regulationDF=regulation.regulation2df(regulation.getRegulation(chrpos[i][0],chrpos[i][1]))

        LOGGER.info("Looking for GWAS hits around the variant")
        gwasDF=gwas.hitsByRegion(chrpos[i][0],chrpos[i][1])
        
        LOGGER.info("Creating VEP dataframe")
        temp=vep.getVepDF(variant_data["rs"],mappings,build=build)
        vepDFtr=temp["transcript"]
        vepDFreg=temp["regulatory"]
        
        LOGGER.info("Creating population dataframe")
        populationDF=variant.population2df(variant_data["population_data"])
        tmp_name=utils.df2svg(populationDF,current_variant_id)
        populationFname=None
        if not tmp_name is None:
            populationFname="./"+os.path.relpath(tmp_name,config.OUTPUT_DIR)
        
        LOGGER.info("Creating PubMed dataframe")
        pubmedDF=pubmed.getPubmedDF(current_variant_id,variant_data["synonyms"])
        
        LOGGER.info("Retrieving genes around the variant")
        gene_list=gene.getGeneList(chrpos[i][0],chrpos[i][1],build=build)
        if gene_list:
            all_genes.extend(gene_list)
            LOGGER.info("Got %d gene(s)" %(len(gene_list)))

        LOGGER.info("Creating gene dataframe")
        geneDF=gene.geneList2df(gene_list)
        
        LOGGER.info("Creating gnomAD dataframe")
        gnomadDF=gnomad.getPopulationAF(current_variant_id)
        tmp_name=utils.df2svg(gnomadDF,current_variant_id)
        gnomadFname=None
        if not tmp_name is None:
            gnomadFname="./"+os.path.relpath(tmp_name,config.OUTPUT_DIR)
            
        LOGGER.info("Creating GTEx dataframe")
        GTEx_genesDF=gtex.getGTExDF(mappings)    
        LOGGER.info("Found %d eQTL(s)\n" % len(GTEx_genesDF))

        if out_html:
            if len(variantDF)>0:
                D["variant_table%d" %i]=variantDF.to_html(index=True,classes='utf8',table_id="common")
            if len(regulationDF)>0:
                D["regulation_table%d" %i]=regulationDF.to_html(index=False,classes='utf8',table_id="common")
            if len(gwasDF)>0:
                D["gwas_table%d" %i]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
            if len(vepDFtr)>0:
                D["vep_table%d" %i]=vepDFtr.to_html(index=False,classes='utf8',table_id="common")
            if len(vepDFreg)>0:
                D["vepreg_table%d" %i]=vepDFreg.to_html(index=False,classes='utf8',table_id="common")
            if not populationFname is None:
                D["population_table%d" %i]="<img src=\"%s\">" % populationFname
            if not gnomadFname is None:
                D["gnomad_table%d" %i]="<img src=\"%s\">" % gnomadFname
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
            D1["population_table%d" %i]=populationDF.to_json(orient="records")
            D1["pubmed_table%d" %i]=pubmedDF.to_json(orient="records")
            D1["phenotype_table%d" %i]=phenotypeDF.to_json(orient="records")
            D1["gene_table%d" %i]=geneDF.to_json(orient="records")
            D1["gtex_genes_table%d" %i]=GTEx_genesDF.to_json(orient="records")

#----------------------------------------------------------------------------------------------------------------------
            
    # Gene information for all genes near the variant

    gene_names=list()
    LOGGER.info("Total genes: %d" % len(all_genes))
    for i in range(0,len(all_genes)):
        gene_ID=all_genes[i]["ID"]
        gene_names.append(gene_ID)
        LOGGER.info("Current gene: %s",gene_ID)
        LOGGER.info("Retreiving general information")
        info=gene.getGeneInfo(gene_ID,build=build)
        infoDF=gene.geneInfo2df(info)
        
        LOGGER.info("Retreiveing cross-references")
        xrefs=gene.getGeneXrefs(gene_ID)
        LOGGER.info("Creating GO terms dataframe")
        goDF=gene.goterms2df(xrefs)
        if xrefs is None:
            up=None
        else:
            LOGGER.info("Retreiving UniProt data")
            if len(xrefs["UniProtKB/Swiss-Prot"])>0:
                up=uniprot.getUniprotData(xrefs["UniProtKB/Swiss-Prot"][0][0])
            else:
                up=None
        LOGGER.info("Creating UniProt dataframe")
        uniprotDF=uniprot.uniprot2df(up)
        
        LOGGER.info("Retreiving GWAS data")
        gwasDF=gwas.hitsByGene(gene_ID)
        
        LOGGER.info("Retreiving GTEx data")
        gtexDF=gtex.gtex2df(gtex.parseGTEx(info["chromosome"],info["start"],info["end"],gene_ID))
        
        LOGGER.info("Retreiving mouse data")
        mouseDF=mouse.getMousePhenotypes(gene_ID)
        
        LOGGER.info("Creating GXA dataframe")
        gxaDF=gxa.getGxaDFLocal(gene_ID)
        # relative to the output dir
        tmp_name=gxa.df2svg(gxaDF)
        gxaFname=None
        if not tmp_name is None:
            gxaFname="./"+os.path.relpath(tmp_name,config.OUTPUT_DIR)
        LOGGER.info("")

        if out_html:
            if len(infoDF)>0:
                D["gene_table_%s" % gene_ID]=infoDF.to_html(index=False,classes='utf8',table_id="common")
            if len(goDF):
                D["go_table_%s" % gene_ID]=goDF.to_html(index=False,classes='utf8',table_id="common")
            if len(uniprotDF)>0:
                D["uniprot_table_%s" % gene_ID]=uniprotDF.to_html(index=False,classes='utf8',table_id="common")
            if len(gwasDF)>0:
                D["gwas_table_%s" % gene_ID]=gwasDF.to_html(index=False,classes='utf8',table_id="common",render_links=True,escape=False)
            if not gxaFname is None:
                D["gxa_heatmap_%s" % gene_ID]="<img src=\"%s\">" % gxaFname
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
                
#--------------------------------------------- CREATING OUTPUT ---------------------------------------------------------

    if out_html:
        outfile=config.OUTPUT_DIR+"/%s.html" % current_variant_id
        LOGGER.info("Saving HTML output to %s\n" % outfile)
        template_fname=config.OUTPUT_DIR+"/template.html"
        utils.generateTemplate(mapping_names,gene_names,template_fname)
        f = open(outfile,"w")
        f.write(utils.generateHTML(template_fname,D))
        f.close()
        if os.path.isfile(template_fname):
            os.remove(template_fname)
    else:
        outfile=config.OUTPUT_DIR+"/%s.json.gz" %current_variant_id
        LOGGER.info("Saving JSON output to %s\n" % outfile)
        with gzip.GzipFile(outfile,"w") as fout:
            fout.write(json.dumps(D1).encode('utf-8'))

#----------------------------------------------------------------------------------------------------------------------

if os.path.isfile(config.OUTPUT_DIR+"/expression_atlas-homo_sapiens.tsv"):
    os.remove(config.OUTPUT_DIR+"/expression_atlas-homo_sapiens.tsv")

