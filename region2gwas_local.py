#!/usr/bin/python3

import sys
import re
import argparse
import logging
import json
import pandas as pd

import locale
locale.setlocale(locale.LC_CTYPE,"en_US.UTF-8")

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="Get GWAS catalog associations for a given region")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level",required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--region','-r',action="store",help="Region: chr:start-end",required=False,default=None)
    parser.add_argument('--gwas','-g',action="store",help="GWAS catalog file",required=True,default=None)
    parser.add_argument('--json','-j', action="store_true",help="Output JSON",required=False)

    try:
        args=parser.parse_args()
    except ValueError as e:
        print(str(e))
        parser.print_help()
        sys.exit(1)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    gwas=args.gwas
    region=args.region
    out_json=args.json

    LOGGER=logging.getLogger("region2gwas_local")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    if sys.stdin.isatty() and region is None:
        parser.print_help()
        sys.exit(1)
        
    m=re.match("(\w+):(\d+)-(\d+)",region)
    if m:
        chrom=m.group(1)
        start=int(m.group(2))
        end=int(m.group(3))
    else:
        LOGGER.error("Malformed region")
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    try:
        T=pd.read_table(gwas,sep="\t",dtype={"CHR_ID":str,"CHR_POS":int,"SNPS":str,"P-VALUE":str,"DISEASE/TRAIT":str,"PUBMEDID":str},encoding="utf-8")
    except Exception as e:
        LOGGER.error(str(e))
        sys.exit(1)

    T.fillna(value="NA",inplace=True)
    df=pd.DataFrame(T[(T["CHR_ID"]==chrom) & (T["CHR_POS"]>start) & (T["CHR_POS"]<end)])
    # df["P-VALUE"]=df["P-VALUE"].astype(float)
    df.sort_values(by="P-VALUE",key=lambda x:x.astype(float),inplace=True)
    if out_json:   
        df[["SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID"]].to_json(sys.stdout,orient="records")
    else:
        df[["SNPS","P-VALUE","DISEASE/TRAIT","PUBMEDID"]].to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",encoding="utf-8")
            
if __name__=="__main__":
    main()
