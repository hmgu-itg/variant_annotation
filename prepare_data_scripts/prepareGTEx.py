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
parser.add_argument('--temp','-t',default="/tmp/",action="store_true",help="Temp dir")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input','-i', action="store",help="input GTEx.tar",required=True)
requiredArgs.add_argument('--gencode','-g', action="store",help="input Gencode gtf.gz file",required=True)
requiredArgs.add_argument('--lookup','-l', action="store",help="input GTEx lookup txt.gz file",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

version="8"
tdir="/tmp/"
if args.version:
    version=args.version
if args.temp:
    tdir=args.temp

if not tdir.endswith("/"):
    tdir+="/"

fname=args.input
gname=args.gencode
lname=args.lookup

# ======================================== BUILDING ID MAPPING ======================================================

id_mapping=dict() # ID --> rsID

if version=="8":
    ldf=pd.read_table(lname,header=0,compression="gzip")
    for index, row in ldf.iterrows():
        var=row["variant_id"]
        if var.startswith("chr"):
            var=var[3:]
        if var.endswith("_b38"):
            var=var[:-4]
        rs=row["rs_id_dbSNP151_GRCh38p7"]
        if var in id_mapping:
            print("WARNING: %s occurs multiple times in %s" % (var,lname),file=sys.stderr)
            continue
        else:
            id_mapping[var]=rs

# ======================================= BUILDING GENE MAPPING =====================================================

gene_map=dict() # gene --> chr,start,end

if version=="8":
    gdf=pd.read_table(gname,header=None,comment="#",compression="gzip")
    gdf=gdf[gdf[2]=="gene"]
    for index, row in gdf.iterrows():
        c=row[0]
        if c.startswith("chr"):
            c=c[3:]
        s=row[3]
        e=row[4]
        f=row[8]
        m=re.search("gene_id\s+\"(ENSG\d+)(\.\d+)?(\_PAR\_Y)?\"",f)
        if m.group(3):
            print("WARNING: %s skipped (PAR Y)." % (m.group(0)),file=sys.stderr)
            continue
        if m:
            gene_id=m.group(1)
            if gene_id in gene_map:
                print("WARNING: %s occurs multiple times in %s" % (gene_id,gname),file=sys.stderr)
                continue
            else:
                gene_map[gene_id]={"chr":c,"start":int(s),"end":int(e)}
        else:
            print("ERROR: no gene ID in row %d" % index,file=sys.stderr)
            print(row)
            continue

# ======================================== BUILDING VAR MAPPING ======================================================

    gene2var=dict()
    var_map=dict()

    infile=tf.open(fname)
    files=infile.getmembers()
    for f in files:
        m=re.search("/([^.]+)\.v8\.signif_variant_gene_pairs\.txt\.gz",f.name)
        if m:
            tissue=m.group(1)
            infile.extract(f,tdir)
            T=pd.read_table(tdir+f.name,compression="gzip",header=0)
            for index, row in T.iterrows():
                var=row["variant_id"]
                if var.endswith("_b38"):
                    var=var[:-4]
                if var.startswith("chr"):
                    var=var[3:]
                if var not in var_map:
                    var_map[var]={tissue:[]}
                elif tissue not in var_map[var]:
                    var_map[var][tissue]=[]

                gene=row["gene_id"]
                pval=row["pval_nominal"]
                beta=row["slope"]
                SE=row["slope_se"]
                dist=row["tss_distance"]
                m1=re.search("^(ENSG\d+)",gene)
                if m1:
                    gene=m1.group(1)
                    if tissue not in gene2var:
                        gene2var[tissue]={gene:[]}
                    elif gene not in gene2var[tissue]:
                        gene2var[tissue][gene]=[]

                    if next((z for z in var_map[var][tissue] if z["gene"]==gene),None):
                        print("WARNING: in file",f.name,"variant",var,"has multiple  associations with",gene,sep=" ",file=sys.stderr)
                    else:
                        var_map[var][tissue].append({"gene":gene,"p":pval,"beta":beta,"SE":SE,"dist":dist})
                        if var in gene2var[tissue][gene]:
                            print("WARNING: in file %s variant %s has multiple associations with %s" %(f.name,var,gene),file=sys.stderr)
                        else:
                            gene2var[tissue][gene].append(var)

            print(tissue,file=sys.stderr)
            print(f.name,file=sys.stderr)
            sys.stderr.flush()

    infile.close()

# =============================================== OUTPUT ===============================================================
    
    for tissue in gene2var:
        for gene in gene2var[tissue]:
            if gene in gene_map:
                for var in gene2var[tissue][gene]:
                    rs="."
                    if var in id_mapping:
                        rs=id_mapping[var]
                    else:
                        print("WARNING: no rs for %s" %var,file=sys.stderr)
                    L=var_map[var][tissue]
                    x=next((z for z in L if z["gene"]==gene),None)
                    if x:
                        p=x["p"]
                        beta=x["beta"]
                        SE=x["SE"]
                        dist=x["dist"]
                        print("%s\t%d\t%d\t%s\t%s/%s:%s:%s:%s:%s:%s" % (gene_map[gene]["chr"],gene_map[gene]["start"]-1,gene_map[gene]["end"],gene,var,rs,tissue,p,beta,SE,dist),file=sys.stdout)
                    else:
                        print("ERROR: could not find (%s %s %s) association" %(var, tissue, gene),file=sys.stderr)
            else:
                print("ERROR: gene coordinates for",gene,"were not found",sep=" ",file=sys.stderr)

    for var in var_map:
        rs="."
        if var in id_mapping:
            rs=id_mapping[var]
        else:
            print("WARNING: no rs for %s" %var,file=sys.stderr)
        m=re.search("([^\W_]+)_(\d+)_([ACGT]+)_([ACGT]+)$",var)
        if not m:
            print("ERROR: malformed variant id %s" %(var),file=sys.stderr)
            continue
        c=m.group(1)
        pos=int(m.group(2))
        ref=m.group(3)
        alt=m.group(4)
        start=pos-1
        end=pos+max(len(ref),len(alt))-1
        for t in var_map[var]:
            L=var_map[var][t]
            for x in L:
                print("%s\t%d\t%d\t%s/%s\t%s:%s:%s:%s:%s:%s" %(c,start,end,var,rs,x["gene"],t,x["p"],x["beta"],x["SE"],x["dist"]),file=sys.stdout)
