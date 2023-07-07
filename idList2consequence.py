#!/usr/bin/python3

import sys
import os
import argparse
import json
import logging

from varannot import utils
from varannot import query
from varannot import config
from varannot import variant

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    build="38"
    verbosity=logging.INFO

    parser = argparse.ArgumentParser(description="Get functional consequences for a list of variant IDs")
    parser.add_argument('--build','-b',action="store",help="Genome build: default: 38", default="38",required=False)
    parser.add_argument('--input','-i',action="store",help="Input ID list",required=True)
    parser.add_argument('--name','-n',action="store_true",help="Use gene name instead of ENSEMBL ID",required=False)
    parser.add_argument('--max','-m',action="store_true",help="Only report the most severe consequence",required=False)
    parser.add_argument("--verbose", "-v",help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")

    try:
        args=parser.parse_args()
    except:
       sys.exit(0)

    if args.build!=None:
        build=args.build

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    infile=args.input
    use_name=args.name
    use_max=args.max

    LOGGER=logging.getLogger("idList2consequence")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.utils").addHandler(ch)
    logging.getLogger("varannot.utils").setLevel(verbosity)
    logging.getLogger("varannot.variant").addHandler(ch)
    logging.getLogger("varannot.variant").setLevel(verbosity)
    logging.getLogger("varannot.query").addHandler(ch)
    logging.getLogger("varannot.query").setLevel(verbosity)

    #---------------------------------------------------------------------------------------------------------------------------
    
    IDs=list()
    with open(infile) as f:
        IDs=f.read().splitlines()

    chunks=list()
    i=1
    for L in utils.chunks(IDs,config.VEP_POST_MAX):
        LOGGER.info("Current chunk: %d" % i)
        i+=1
        if use_name:
            r=variant.addConsequencesToIDList(L,build=build,most_severe_only=use_max,gene_key="gene_symbol")
        else:
            r=variant.addConsequencesToIDList(L,build=build,most_severe_only=use_max,gene_key="gene_id")
        if r is None:
            continue
        else:
            chunks.append(r)
    for ID in IDs:
        found=False
        for c in chunks:
            if ID in c:
                found=True
                for g in c[ID]:
                    print("%s\t%s\t%s" %(ID,g,c[ID][g]))
                break
        if not found:
            print("%s\tNA\tNA" %(ID))

if __name__=="__main__":
    main()
