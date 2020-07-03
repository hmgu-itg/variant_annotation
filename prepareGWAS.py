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
for index, row in T[["SNPS","PUBMEDID","CHR_ID","CHR_POS","DISEASE/TRAIT","P-VALUE"]].iterrows():
    a=row["CHR_ID"].split(";")
    b=row["CHR_POS"].split(";")
    c=row["SNPS"].split(";")
    if len(a)==1 and len(b)==1 and len(c)==1:
        print("%s\t%s\t%s\t%s\t%s\t%s" %(row["CHR_ID"],row["CHR_POS"],row["SNPS"],row["P-VALUE"],row["DISEASE/TRAIT"],row["PUBMEDID"]))
    # for f in files:
    #     m=re.search("/([^.]+)\.v8\.egenes.txt.gz",f.name)
    #     if m:
    #         tissue=m.group(1)
    #         infile.extract(f,"/tmp")
    #         T=pd.read_table("/tmp/"+f.name,compression="gzip",header=0)
    #         for index, row in T.iterrows():
    #             id_mapping[row["variant_id"]]=row["rs_id_dbSNP151_GRCh38p7"]
    #             gene=row["gene_id"]
    #             m1=re.search("^(ENSG\d+)",gene)
    #             if m1:
    #                 gene=m1.group(1)
    #                 c1=row["gene_chr"]
    #                 if c1.startswith("chr"):
    #                     c1=c1[3:]
    #                 s1=int(row["gene_start"])
    #                 e1=int(row["gene_end"])
    #                 if gene in gene_map:
    #                     c=gene_map[gene]["chr"]
    #                     s=int(gene_map[gene]["start"])
    #                     e=int(gene_map[gene]["end"])
    #                     if c!=c1 or s!=s1 or e!=e1:
    #                         print("WARNING: in",f.name,": gene",gene,"has multiple coordinates:",sep=" ",file=sys.stderr)
    #                         print(c1,str(s1),str(e1),sep=" ",file=sys.stderr)
    #                         print(c,str(s),str(e),sep=" ",file=sys.stderr)
    #                         print("",file=sys.stderr)
    #                 else:
    #                     gene_map[gene]={"chr":c1,"start":s1,"end":e1}

    #         print(tissue,file=sys.stderr)
    #         print(f.name,file=sys.stderr)

