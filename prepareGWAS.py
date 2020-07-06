#!/usr/bin/python3

import argparse
import sys
import re
import pandas as pd
import functions

parser = argparse.ArgumentParser(description="Clean up GWAS associations file")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input','-i', action="store",help="input GWAS.tsv.gz",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

fname=args.input

T=pd.read_table(fname,compression="gzip",header=0,sep="\t",keep_default_na=False,dtype=str)
print("CHR_ID","CHR_POS","SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID",sep="\t")
for index, row in T[["SNPS","PUBMEDID","CHR_ID","CHR_POS","DISEASE/TRAIT","P-VALUE"]].iterrows():
    a=row["CHR_ID"].split(";")
    b=row["CHR_POS"].split(";")
    c=row["SNPS"].split(";")

    if len(row["CHR_ID"])==0 or len(row["CHR_POS"])==0:
        if len(c)==1 and len(row["SNPS"])!=0:
            m=re.match("(rs\d+)",row["SNPS"],re.I)
            if m:
                rs=m.group(1)
                L=functions.rs2position(rs)
                if L:
                    for z in L:
                        print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
            else:
                m=re.match("^snp(\d+)-(\d+)",row["SNPS"],re.I)
                if m:
                    c=m.group(1)
                    p=m.group(2)
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                else:
                    m=re.match("^(chr)?(\d+):\s*(\d+)",row["SNPS"],re.I)
                    if m:
                        c=m.group(2)
                        p=m.group(3)
                        print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    else:
                        m=re.match("^chr(\d+)[._](\d+)",row["SNPS"],re.I)
                        if m:
                            c=m.group(1)
                            p=m.group(2)
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                        else:
                            m=re.match("^(\d{1,2})-(\d+)",row["SNPS"],re.I)
                            if m:
                                c=m.group(1)
                                p=m.group(2)
                                print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                            else:
                                m=re.match("^([XY]):(\d+)",row["SNPS"],re.I)
                                if m:
                                    c=m.group(1)
                                    p=m.group(2)
                                    print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                                else:
                                    m=re.match("^del-(\d+):(\d+)",row["SNPS"],re.I)
                                    if m:
                                        c=m.group(1)
                                        p=m.group(2)
                                        print("%s\t%s\t%s\t%s\t%s\t%s" %(c,p,row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
    else:
        if len(a)==1 and len(b)==1 and len(c)==1:
            m1=re.match("^(\w+)\s+x\s+(\w+)",row["CHR_ID"],re.I)
            m2=re.match("^(\d+)\s+x\s+(\d+)",row["CHR_POS"],re.I)
            m3=re.match("^(rs\d+)\s+x\s+(rs\d+)",row["SNPS"],re.I)
            if not m1 and not m3:
                print("%s\t%s\t%s\t%s\t%s\t%s" %(row["CHR_ID"],row["CHR_POS"],row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
            else:
                if m2:
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(m1.group(1),m2.group(1),m3.group(1),row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    print("%s\t%s\t%s\t%s\t%s\t%s" %(m1.group(2),m2.group(2),m3.group(2),row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                else:
                    rs=m3.group(1)
                    L=functions.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
                    rs=m3.group(2)
                    L=functions.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
        elif len(c)!=1:
            for v in c:
                m=re.match("(rs\d+)",v,re.I)
                if m:
                    rs=m.group(1)
                    L=functions.rs2position(rs)
                    if L:
                        for z in L:
                            print("%s\t%s\t%s\t%s\t%s\t%s" %(z["chr"],z["pos"],rs,row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))

