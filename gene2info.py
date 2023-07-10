#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd

from varannot import gene

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser = argparse.ArgumentParser(description="General gene info")
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
    ID=args.id
    out_json=args.json

    LOGGER=logging.getLogger("gene2info")
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

    data=gene.getGeneInfo(ID,build=build)
    if data is None:
        sys.exit(1)
    df=pd.json_normalize([data])
    if out_json:
        df.fillna(value="NA",inplace=True)
        df[["id","name","chromosome","start","end","description","type","strand"]].to_json(sys.stdout,orient="records")
    else:
        df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["id","name","chromosome","start","end","description","type","strand"])

if __name__=="__main__":
    main()
