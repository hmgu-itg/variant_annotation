#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd
from functools import partial

from varannot import gene

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="List gene transcripts")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--json','-j', action="store_true",help="Output JSON",required=False)
    parser.add_argument('--id','-i', action="store",help="Gene ID",required=False,default=None)
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38",required=False,default="38")

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    build=args.build
    out_json=args.json
    ID=args.id

    LOGGER=logging.getLogger("gene2transcripts")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.gene").addHandler(ch)
    logging.getLogger("varannot.gene").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    data=gene.getGeneTranscripts(ID,build=build)
    if data is None:
        sys.exit(1)
    df=pd.json_normalize(data)
    df.rename(columns={"Parent":"Gene ID","seq_region_name":"chr"},inplace=True)
    df["translation"],df["aa"],df["uniprot"]=zip(*df["id"].map(partial(gene.getTranslationInfo,build=build)))
    if out_json:
        df.fillna(value="NA",inplace=True)
        df=df.astype(str,copy=True)
        df["biotype"]=df["biotype"].apply(lambda x:x.replace("_"," "))
        # df[["Gene ID","id","chr","start","end","biotype","translation","aa","uniprot"]].to_json(sys.stdout,orient="records")
        df[["id","biotype","translation","aa","uniprot"]].to_json(sys.stdout,orient="records")
    else:
        df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["Gene ID","id","chr","start","end","biotype","translation","aa","uniprot"])

if __name__=="__main__":
    main()
