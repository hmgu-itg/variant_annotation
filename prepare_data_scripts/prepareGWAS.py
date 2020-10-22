#!/usr/bin/python3

import logging
import argparse
import sys
import re
import pandas as pd
import os

from varannot import variant
from varannot import utils
from varannot import config

sys.stdout=open(sys.stdout.fileno(),mode='w',encoding='utf8',buffering=1)

verbosity=logging.INFO

parser = argparse.ArgumentParser(description="Clean up GWAS associations file")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input','-i', action="store",help="input gwas_full.tsv.gz",required=True)
requiredArgs.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(0)

fname=args.input

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("prepareGWAS")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.variant").addHandler(ch)
logging.getLogger("varannot.variant").setLevel(verbosity)
logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(verbosity)

#=========================================================================================================================================================

L=list() # list of dicts (from input): chr:pos:id:pval:trait:pmid
T=pd.read_table(fname,header=0,sep="\t",keep_default_na=False,dtype=str,encoding="utf-8",compression="gzip")

LOGGER.debug("Start analyzing input")
count=0

for index, row in T[["SNPS","PUBMEDID","CHR_ID","CHR_POS","DISEASE/TRAIT","P-VALUE"]].iterrows():
    chrID=row["CHR_ID"]
    chrPOS=row["CHR_POS"]
    SNPs=row["SNPS"]
    trait=row["DISEASE/TRAIT"]
    pmid=row["PUBMEDID"]
    pval=row["P-VALUE"]

    count+=1
    if count%10000==0:
        LOGGER.debug("Record %d (%d)" % (count,len(T)))
    
    # case 1: there is only one chr:pos, but several SNPs, because a haplotype is reported
    m1=re.search("^\d+$",chrPOS)
    m2=re.search("[;,x]",SNPs)
    if m1 and m2:
        for x in re.split("\s*[;,x]\s*",SNPs):
            d={"chr":chrID,"pos":chrPOS,"id":x,"trait":trait,"pval":pval,"pmid":pmid}
            if not d in L:
                L.append(d)
        continue

    # case 2: only one field in CHR, no separators in SNPS
    m1=re.search("^[XY\d]+$",chrID,re.I)
    if m1 and not m2:
        d={"chr":chrID,"pos":chrPOS,"id":SNPs,"trait":trait,"pval":pval,"pmid":pmid}
        if not d in L:
            L.append(d)
        continue

    # case 3: CHR has several fields separated by 'x' or ';'
    m1=re.search("\s*[x;]\s*",chrID)
    if m1:
        L1=re.split("\s*[x;]\s*",chrID)
        L2=re.split("\s*[x;]\s*",chrPOS)
        L3=re.split("\s*[x;]\s*",SNPs)
        # TODO: check lengths
        for (c,p,x) in zip(L1,L2,L3):
            d={"chr":c,"pos":p,"id":x,"trait":trait,"pval":pval,"pmid":pmid}
            if not d in L:
                L.append(d)
        continue

    # case 4: both CHR and POS are empty
    m1=re.search("^\s*$",chrID)
    m2=re.search("^\s*$",chrPOS)
    if m1 and m2:
        m3=re.search("(rs\d+)",SNPs)
        if m3:
            d={"chr":"NA","pos":"NA","id":m3.group(1),"trait":trait,"pval":pval,"pmid":pmid}
            if not d in L:
                L.append(d)
        continue
            
        m3=re.search("\*",SNPs)
        if not m3:
            m4=re.search("([XY\d]+)[:_.-]\s*(\d+)",SNPs,re.I)
            if m4:
                d={"chr":m4.group(1),"pos":m4.group(2),"id":SNPs,"trait":trait,"pval":pval,"pmid":pmid}
                if not d in L:
                    L.append(d)
            continue

LOGGER.debug("Done reading input")
#=========================================================================================================================================================

# liftOver records with chr:pos present; if chr:pos==NA:NA, make a REST query

liftover_in=list()
rest_in=set()
LOGGER.debug("Preparing data")
for x in L:
    if x["chr"]=="NA" and x["pos"]=="NA":
        rest_in.add(x["id"])
    else:
        c=x["chr"]
        if not c.startswith("chr"):
            c="chr"+c
        liftover_in.append({"chr":c,"start":int(x["pos"])-1,"end":int(x["pos"]),"id":x["id"]})
LOGGER.debug("Done")

LOGGER.debug("Start liftOver")
liftover_out=utils.runLiftOver(liftover_in,build="37")
LOGGER.debug("Done")

rest_out=dict()
LOGGER.debug("Start REST (%d records)" % len(rest_in))
count=0
for chunk in utils.chunks(list(rest_in),config.VARIANT_RECODER_POST_MAX):
    res=variant.rsList2position(chunk,build="38",alleles=False)
    for x in res:
        rest_out[x]=res[x]
    count+=1
    LOGGER.debug("REST chunk 1")
LOGGER.debug("Done")

#============================================================= OUTPUT ====================================================================================

print("CHR_ID","CHR_POS","SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID",sep="\t")
for x in L:
    if x["chr"]=="NA" and x["pos"]=="NA":
        if x["id"] in rest_out:
            for a in rest_out[x["id"]]:
                print(a["chr"],a["pos"],x["id"],x["pval"],x["trait"],x["pmid"])
    else:
        z=next((a for a in liftover_out if a["id"]==x["id"]),None)
        if z:
            print(print(z["chr"],z["pos"],x["id"],x["pval"],x["trait"],x["pmid"]))
