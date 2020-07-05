#!/usr/bin/python3

import tarfile as tf
import argparse
import sys
import re
import io
import gzip
import pandas as pd

parser = argparse.ArgumentParser(description="Prepare GTEx BED file")
parser.add_argument('--version','-v',default="8",action="store_true",help="GTEx version")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input','-i', action="store",help="input GTEx.tar",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

version="8"
if args.version:
    version=args.version

fname=args.input

# ======================================== BUILDING ID AND GENE MAPPING ======================================================

id_mapping=dict() # ID --> rsID
gene_map=dict() # gene --> chr,start,end

if version=="8":
    infile=tf.open(fname)
    names=infile.getnames()
    files=infile.getmembers()
    for f in files:
        m=re.search("/([^.]+)\.v8\.egenes.txt.gz",f.name)
        if m:
            tissue=m.group(1)
            infile.extract(f,"/tmp")
            T=pd.read_table("/tmp/"+f.name,compression="gzip",header=0)
            for index, row in T.iterrows():
                id_mapping[row["variant_id"]]=row["rs_id_dbSNP151_GRCh38p7"]
                gene=row["gene_id"]
                m1=re.search("^(ENSG\d+)",gene)
                if m1:
                    gene=m1.group(1)
                    c1=row["gene_chr"]
                    if c1.startswith("chr"):
                        c1=c1[3:]
                    s1=int(row["gene_start"])
                    e1=int(row["gene_end"])
                    if gene in gene_map:
                        c=gene_map[gene]["chr"]
                        s=int(gene_map[gene]["start"])
                        e=int(gene_map[gene]["end"])
                        if c!=c1 or s!=s1 or e!=e1:
                            print("WARNING: in",f.name,": gene",gene,"has multiple coordinates:",sep=" ",file=sys.stderr)
                            print(c1,str(s1),str(e1),sep=" ",file=sys.stderr)
                            print(c,str(s),str(e),sep=" ",file=sys.stderr)
                            print("",file=sys.stderr)
                    else:
                        gene_map[gene]={"chr":c1,"start":s1,"end":e1}

            print(tissue,file=sys.stderr)
            print(f.name,file=sys.stderr)

# ======================================== BUILDING VAR MAPPING ======================================================

    gene2var=dict()
    var_map=dict()

    for f in files:
        m=re.search("/([^.]+)\.v8\.signif_variant_gene_pairs\.txt\.gz",f.name)
        if m:
            tissue=m.group(1)
            infile.extract(f,"/tmp")
            T=pd.read_table("/tmp/"+f.name,compression="gzip",header=0)
            for index, row in T.iterrows():
                var=row["variant_id"]
                if var not in var_map:
                    var_map[var]={tissue:[]}
                elif tissue not in var_map[var]:
                    var_map[var][tissue]=[]
                gene=row["gene_id"]
                pval=row["pval_nominal"]
                beta=row["slope"]
                SE=row["slope_se"]
                m1=re.search("^(ENSG\d+)",gene)
                if m1:
                    gene=m1.group(1)
                    if next((z for z in var_map[var][tissue] if z["gene"]==gene),None):
                        print("WARNING: in file",f.name,"variant",var,"has multiple  associations with",gene,sep=" ",file=sys.stderr)
                    else:
                        var_map[var][tissue].append({"gene":gene,"p":pval,"beta":beta,"SE":SE})
                        if gene in gene2var:
                            gene2var[gene].append(var)
                        else:
                            gene2var[gene]=[var]

            print(tissue,file=sys.stderr)
            print(f.name,file=sys.stderr)

    infile.close()

# =============================================== OUTPUT ===============================================================
    
    for gene in gene2var:
        if gene in gene_map:
            for var in gene2var[gene]:
                a=var_map[var]
                for t in a:
                    x=next((z for z in a[t] if z["gene"]==gene),None)
                    if x:
                        p=x["p"]
                        beta=x["beta"]
                        SE=x["SE"]
                        print("%s\t%d\t%d\t%s\t%s:%s:%s:%s:%s" % (gene_map[gene]["chr"],gene_map[gene]["start"]-1,gene_map[gene]["end"],gene,var,t,p,beta,SE),file=sys.stdout)
        else:
            print("ERROR: gene coordinates for",gene,"were not found",sep=" ",file=sys.stderr)

    for var in var_map:
        if var in id_mapping:
            id2=id_mapping[var]
        else:
            id2="."
        m=re.search("^(chr)?([^\W_]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38$",var)
        if not m:
            print("ERROR: malformed variant id %s" %(var),file=sys.stderr)
            continue
        c=m.group(2)
        pos=int(m.group(3))
        ref=m.group(4)
        alt=m.group(5)
        start=pos-1
        end=pos+max(len(ref),len(alt))-1
        for t in var_map[var]:
            L=var_map[var][t]
            for x in L:
                print("%s\t%d\t%d\t%s\t%s:%s:%s:%s:%s" %(c,start,end,var,x["gene"],t,x["p"],x["beta"],x["SE"]),file=sys.stdout)
