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

# --------------------------- GWAS Catalog has b38 coordinates ---------------------------

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

T=pd.read_table(fname,header=0,sep="\t",keep_default_na=False,dtype=str,encoding="utf-8",compression="gzip")
total=len(T)
id2cp=dict() # ID --> chr:pos, to track ambiguous IDs
ambiguous=set()
S=set()

LOGGER.debug("Start analyzing input")

for index, row in T[["SNPS","PUBMEDID","CHR_ID","CHR_POS","DISEASE/TRAIT","P-VALUE"]].iterrows():
    chrID=row["CHR_ID"]
    chrPOS=row["CHR_POS"]
    SNPs=row["SNPS"]
    trait=row["DISEASE/TRAIT"]
    pmid=row["PUBMEDID"]
    pval=row["P-VALUE"]

    # case 1: there is only one chr:pos, but several SNPs, because a haplotype is reported
    m1=re.search("^\d+$",chrPOS)
    m2=re.search("[;,x]",SNPs)
    if m1 and m2:
        for x in re.split("\s*[;,x]\s*",SNPs):
            S.add(";".join(["NA","NA",x,pval,pmid,trait]))
        continue

    # case 2: only one field in CHR, no separators in SNPS
    m1=re.search("^[XY\d]+$",chrID,re.I)
    if m1 and not m2:
        S.add(";".join([chrID,chrPOS,SNPs,pval,pmid,trait]))
        s=chrID+":"+chrPOS
        if not SNPs in id2cp:
            id2cp[SNPs]=s
        elif s!=id2cp[SNPs]:
            ambiguous.add(SNPs)
        continue

    # case 3: CHR has several fields separated by 'x' or ';'
    m1=re.search("\s*[x;]\s*",chrID)
    if m1:
        L1=re.split("\s*[x;]\s*",chrID)
        L2=re.split("\s*[x;]\s*",chrPOS)
        L3=re.split("\s*[x;]\s*",SNPs)
        
        if len(L1)!=len(L2) or len(L1)!=len(L3):
            LOGGER.warning("Input record for %s %s %s is malformed" % (chrID,chrPOS,SNPs))
            for x in L3:
                S.add(";".join(["NA","NA",x,pval,pmid,trait]))
            continue
        
        for (c,p,x) in zip(L1,L2,L3):
            S.add(";".join([c,p,x,pval,pmid,trait]))
            s=c+":"+p
            if not x in id2cp:
                id2cp[x]=s
            elif s!=id2cp[x]:
                ambiguous.append(x)
        continue

    # case 4: both CHR and POS are empty
    m1=re.search("^\s*$",chrID)
    m2=re.search("^\s*$",chrPOS)
    m3=re.search("(rs\d+)",SNPs)
    if m1 and m2 and m3:
        S.add(";".join(["NA","NA",m3.group(1),pval,pmid,trait]))
            

LOGGER.debug("Done analyzing input")
#=========================================================================================================================================================

id2cp.clear()

# all ambiguous IDs should be queried
LOGGER.debug("Found %d IDs with ambiguous chr:pos" % len(ambiguous))
to_remove=list()
to_replace=list()
for x in S:
    (c,pos,ID,p,pmid,trait)=x.split(";",maxsplit=5)
    if ID in ambiguous:
        to_remove.append(x)
        to_replace.append(";".join(["NA","NA",ID,p,pmid,trait]))

for (r1,r2) in zip(to_remove,to_replace):
    S.remove(r1)
    S.add(r2)
    
#=========================================================================================================================================================

# if chr:pos==NA:NA, make a REST query

D=dict()
rest_in=set()
LOGGER.debug("Preparing data")
for x in S:
    (c,pos,ID,p,pmid,trait)=x.split(";",maxsplit=5)
    if c=="NA":
        rest_in.add(ID)
    else:
        # TODO: check uniqueness of ID --> chr:pos
        D[ID]=":".join([c,pos])
LOGGER.debug("Done preparing data")

rest_out=dict()
LOGGER.debug("Start REST (%d records)" % len(rest_in))
count=1
for chunk in utils.chunks(list(rest_in),config.VARIANT_RECODER_POST_MAX):
    res=variant.rsList2position(chunk,build="38",alleles=False)
    if res:
        LOGGER.debug("Done REST chunk %d" % count)
        for x in res:
            rest_out[x]=res[x]
    else:
        LOGGER.error("REST query for chunk %d failed" % count)
    count+=1
    
LOGGER.debug("Done, REST output: %d records" % len(rest_out))

#============================================================= OUTPUT ====================================================================================

LOGGER.debug("Start output")
print("CHR_ID","CHR_POS","SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID",sep="\t")
failed_rest=set()
for x in S:
    (c,pos,ID,pval,pmid,trait)=x.split(";",maxsplit=5)
    if c=="NA":
        if ID in rest_out:
            for a in rest_out[ID]:
                print("%s\t%s\t%s\t%s\t%s\t%s" % (a["chr"],a["pos"],ID,pval,trait,pmid))
        else:
            failed_rest.add(ID)
    else:
        (c,pos)=D[ID].split(":")
        print("%s\t%s\t%s\t%s\t%s\t%s" % (c,pos,ID,pval,trait,pmid))
            
LOGGER.debug("Failed in REST query: %d" % len(failed_rest))
LOGGER.debug("Done")
