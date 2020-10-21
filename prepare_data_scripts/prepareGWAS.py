#!/usr/bin/python3

import logging
import argparse
import sys
import re
import pandas as pd
import os

from varannot import variant
from varannot import utils

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
print("CHR_ID","CHR_POS","SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID",sep="\t")
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
    # if (/([XY\d]+)[:_.-]\s*(\d+)/i && !/\*/)
    if m1 and m2:
        m3=re.search("\*",SNPs)
        if not m3:
            m4=re.search("([XY\d]+)[:_.-]\s*(\d+)",SNPs,re.I)


    if len(row["CHR_ID"])==0 or len(row["CHR_POS"])==0:
        if len(c)==1 and len(row["SNPS"])!=0:
            m=re.search("(rs\d+)",row["SNPS"],re.I)
            if m:
                rs=m.group(1)
                L=variant.rs2position(rs)
                if L:
                    for z in L:
                        print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
            else:
                m=re.search("^snp(\d+)-(\d+)",row["SNPS"],re.I)
                if m:
                    c=m.group(1)
                    p=m.group(2)
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                else:
                    m=re.search("^(chr)?(\d+):\s*(\d+)",row["SNPS"],re.I)
                    if m:
                        c=m.group(2)
                        p=m.group(3)
                        print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    else:
                        m=re.search("^chr(\d+)[._](\d+)",row["SNPS"],re.I)
                        if m:
                            c=m.group(1)
                            p=m.group(2)
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                        else:
                            m=re.search("^(\d{1,2})-(\d+)",row["SNPS"],re.I)
                            if m:
                                c=m.group(1)
                                p=m.group(2)
                                print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                            else:
                                m=re.search("^([XY]):(\d+)",row["SNPS"],re.I)
                                if m:
                                    c=m.group(1)
                                    p=m.group(2)
                                    print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                                else:
                                    m=re.search("^del-(\d+):(\d+)",row["SNPS"],re.I)
                                    if m:
                                        c=m.group(1)
                                        p=m.group(2)
                                        print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
    else:
        if len(a)==1 and len(b)==1 and len(c)==1:
            m1=re.search("^(\w+)\s+x\s+(\w+)",row["CHR_ID"],re.I)
            m2=re.search("^(\d+)\s+x\s+(\d+)",row["CHR_POS"],re.I)
            m3=re.search("^(rs\d+)\s+x\s+(rs\d+)",row["SNPS"],re.I)
            if not m1 and not m3:
                print("%s\t%s\t%s\t%s\t%s\t%s" %(row["CHR_ID"],row["CHR_POS"],row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
            else:
                if m2:
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(m1.group(1),m2.group(1),m3.group(1),row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(m1.group(2),m2.group(2),m3.group(2),row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                else:
                    rs=m3.group(1)
                    L=variant.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    rs=m3.group(2)
                    L=variant.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
        elif len(c)!=1:
            for v in c:
                m=re.search("(rs\d+)",v,re.I)
                if m:
                    rs=m.group(1)
                    L=variant.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))

