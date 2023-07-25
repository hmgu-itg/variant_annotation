#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging

from varannot import variant
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    build="38"
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="Get rs ID for a list of variant ID")
    parser.add_argument('--build','-b',action="store",help="Genome build: default: 38",default="38",required=False)
    parser.add_argument('--input','-i',action="store",help="Input ID list",required=True)
    parser.add_argument('--rs','-r',action="store_true",help="Only output rs IDs",required=False,default=False)
    parser.add_argument('--unique','-u',action="store_true",help="Only report one ID",required=False,default=False)
    parser.add_argument("--verbose", "-v",help="Optional: verbosity level",required=False,choices=("debug","info","warning","error"),default="info")

    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    try:
        args=parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    build=args.build
    infile=args.input
    rs_only=args.rs
    output_unique=args.unique

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    LOGGER=logging.getLogger("id2rs_list")
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

    #---------------------------------------------------------------------------------------------------------------------------

    IDs=list()
    with open(infile) as f:
        IDs=f.read().splitlines()

    R=variant.id2rs_list(IDs,build)
    if rs_only:
        S1=set()
        for v in R:
            S2=set()
            S=R[v]
            for x in S:
                if x.startswith("rs"):
                    S2.add(x)
            if len(S2)==0:
                S1.add(v)
        for v in S1:
            del R[v]
    else:
        S1=set()
        for v in R:
            if len(R[v])==0:
                S1.add(v)
        for v in S1:
            del R[v]
            
    for v in IDs:
        if v in R:
            if output_unique:
                print("%s\t%s" % (v,list(R[v])[0]))
            else:
                for rs in R[v]:
                    print("%s\t%s" % (v,rs))
        else:
            print("%s\t%s" % (v,v))

if __name__=="__main__":
    main()
