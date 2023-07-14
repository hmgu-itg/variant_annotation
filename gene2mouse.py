#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd

from varannot import mouse

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="MGI mouse phenotypes")
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

    LOGGER=logging.getLogger("gene2mouse")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.mouse").addHandler(ch)
    logging.getLogger("varannot.mouse").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    data=mouse.getMousePhenotypes(ID,build=build)
    if data is None:
        sys.exit(1)
    if out_json:
        data.to_json(sys.stdout,orient="records")
    else:
        data.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA")

if __name__=="__main__":
    main()
