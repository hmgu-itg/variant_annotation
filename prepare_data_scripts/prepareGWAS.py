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

T=pd.read_table(fname,header=0,sep="\t",keep_default_na=False,dtype=str,encoding="utf-8",compression="gzip")
total=len(T)
id2cp=dict() # ID --> chr:pos, to track ambiguous IDs
ambiguous=set()
S=set()

LOGGER.debug("Start analyzing input")
count=0

for index, row in T[["SNPS","PUBMEDID","CHR_ID","CHR_POS","DISEASE/TRAIT","P-VALUE"]].iterrows():
    chrID=row["CHR_ID"]
    chrPOS=row["CHR_POS"]
    SNPs=row["SNPS"]
    trait=row["DISEASE/TRAIT"]
    pmid=row["PUBMEDID"]
    pval=row["P-VALUE"]

    # count+=1
    # if count%10000==0:
    #     LOGGER.debug("Record\t%d/%d" % (count,total))
    
    # case 1: there is only one chr:pos, but several SNPs, because a haplotype is reported
    m1=re.search("^\d+$",chrPOS)
    m2=re.search("[;,x]",SNPs)
    if m1 and m2:
        for x in re.split("\s*[;,x]\s*",SNPs):
            S.add(":".join([chrID,str(int(chrPOS)-1),chrPOS,x,pval,pmid,trait]))
            s=chrID+":"+chrPOS
            if not x in id2cp:
                id2cp[x]=s
            elif s!=id2cp[x]:
                ambiguous.add(x)
        continue

    # case 2: only one field in CHR, no separators in SNPS
    m1=re.search("^[XY\d]+$",chrID,re.I)
    if m1 and not m2:
        S.add(":".join([chrID,str(int(chrPOS)-1),chrPOS,SNPs,pval,pmid,trait]))
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
                S.add(":".join(["NA","NA","NA",x,pval,pmid,trait]))
            continue
        
        for (c,p,x) in zip(L1,L2,L3):
            S.add(":".join([c,str(int(p)-1),p,x,pval,pmid,trait]))
            s=c+":"+p
            if not x in id2cp:
                id2cp[x]=s
            elif s!=id2cp[x]:
                ambiguous.append(x)
        continue

    # case 4: both CHR and POS are empty
    m1=re.search("^\s*$",chrID)
    m2=re.search("^\s*$",chrPOS)
    if m1 and m2:
        m3=re.search("(rs\d+)",SNPs)
        if m3:
            S.add(":".join(["NA","NA","NA",m3.group(1),pval,pmid,trait]))
        continue
            
        m3=re.search("\*",SNPs)
        if not m3:
            m4=re.search("([XY\d]+)[:_.-]\s*(\d+)",SNPs,re.I)
            if m4:
                S.add(":".join([m4.group(1),str(int(m4.group(2))-1),str(int(m4.group(2))),SNPs,pval,pmid,trait]))
            continue

LOGGER.debug("Done analyzing input")
#=========================================================================================================================================================

id2cp.clear()

# all ambiguous IDs should be queried
LOGGER.debug("Found %d IDs with ambiguous chr:pos" % len(ambiguous))
to_remove=list()
to_replace=list()
for x in S:
    (c,s,e,ID,p,pmid,trait)=x.split(":",maxsplit=6)
    if ID in ambiguous:
        to_remove.append(x)
        to_replace.append(":".join(["NA","NA","NA",ID,p,pmid,trait]))

for (r1,r2) in zip(to_remove,to_replace):
    S.remove(r1)
    S.add(r2)
    
#=========================================================================================================================================================

# liftOver records with chr:pos present; if chr:pos==NA:NA, make a REST query

liftover_in=list()
rest_in=set()
LOGGER.debug("Preparing data")
for x in S:
    (c,s,e,ID,p,pmid,trait)=x.split(":",maxsplit=6)
    if c=="NA":
        rest_in.add(ID)
    else:
        if not c.startswith("chr"):
            c="chr"+c
        liftover_in.append({"chr":c,"start":s,"end":e,"id":ID})
LOGGER.debug("Done preparing data")

# no need to liftOver, coordinates are already b38
# LOGGER.debug("Start liftOver")
# liftover_out=utils.runLiftOver(liftover_in,build="37")
# LOGGER.debug("LiftOver Done")

rest_out=dict()
LOGGER.debug("Start REST (%d records)" % len(rest_in))
count=0
for chunk in utils.chunks(list(rest_in),config.VARIANT_RECODER_POST_MAX):
    res=variant.rsList2position(chunk,build="38",alleles=False)
    for x in res:
        rest_out[x]=res[x]
    count+=1
    LOGGER.debug("Done REST chunk %d" % count)
LOGGER.debug("Done, REST output: %d records" % len(rest_out))

#=========================================================================================================================================================

# convert to a dictionary
LOGGER.debug("Converting to a dictionary")
D=dict()
#for x in liftover_out:
for x in liftover_in:
    # if x["id"] in D:
    #     LOGGER.warning("%s occurs multiple times" % x["id"])
    D[x["id"]]=":".join([x["chr"],x["end"]])
LOGGER.debug("Done")

#============================================================= OUTPUT ====================================================================================

LOGGER.debug("Start output")
print("CHR_ID","CHR_POS","SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID",sep="\t")
failed_rest=set()
#failed_liftover=set()
for x in S:
    (c,s,e,ID,pval,pmid,trait)=x.split(":",maxsplit=6)
    if c=="NA":
        if ID in rest_out:
            for a in rest_out[ID]:
                print("%s\t%s\t%s\t%s\t%s\t%s" % (a["chr"],a["pos"],ID,pval,trait,pmid))
        else:
            failed_rest.add(ID)
    else:
        (c,pos)=D[ID].split(":")
        print("%s\t%s\t%s\t%s\t%s\t%s" % (c[c.startswith("chr") and 3:],pos,ID,pval,trait,pmid))
        # if ID in D:
        #     (c,pos)=D[ID].split(":")
        #     print("%s\t%s\t%s\t%s\t%s\t%s" % (c[c.startswith("chr") and 3:],pos,ID,pval,trait,pmid))
        # else:
        #     failed_liftover.add(ID)
            
LOGGER.debug("Failed in REST query: %d" % len(failed_rest))
#LOGGER.debug("Failed in liftOver: %d" % len(failed_liftover))
LOGGER.debug("Done")
