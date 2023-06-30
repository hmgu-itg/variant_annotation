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
    parser = argparse.ArgumentParser(description="Get gene ID for gene name")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--name','-n', action="store",help="Gene name",required=False,default=None)
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
    name=args.name

    LOGGER=logging.getLogger("gene2id")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.gene").addHandler(ch)
    logging.getLogger("varannot.gene").setLevel(verbosity)

    if sys.stdin.isatty() and name is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    r=gene.getGeneID(name,build=build)
    if r is None:
        sys.exit(1)
    for ID in r:
        print("%s\t%s" %(name,ID))

if __name__=="__main__":
    main()
