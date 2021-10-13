#!/usr/bin/python3

import sys
import os
import argparse
import json
import logging

from varannot import utils
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="Get VEP consequences for rs ID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
parser.add_argument('--id','-i', action="store",help="rs ID",required=True)
parser.add_argument('--gene','-g', action="store",help="Gene name/ID",required=False)
parser.add_argument('--name','-n', action="store_true",help="Use gene name instead of ENSEMBL ID",required=False)
parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")

try:
    args=parser.parse_args()
except:
    parser.print_help()
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

rsID=args.id
gene=args.gene

use_name=args.name

LOGGER=logging.getLogger("rs2consequences")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(verbosity)

#---------------------------------------------------------------------------------------------------------------------------

msqs=list()
g2c=dict()
r=query.restQuery(query.makeVepRSQueryURL(rsID,build=build))
LOGGER.debug(json.dumps(r,indent=4,sort_keys=True))

if r:
    for x in r:
        msqs.append(x["most_severe_consequence"])
        if "transcript_consequences" in x:
            for t in x["transcript_consequences"]:
                if use_name:
                    g2c.setdefault(t["gene_symbol"],[]).extend(t["consequence_terms"])
                else:
                    g2c.setdefault(t["gene_id"],[]).extend(t["consequence_terms"])
if gene in g2c:
    msq=utils.getMostSevereConsequence(g2c[gene])
    print("%s\t%s\t%s" %(rsID,gene,msq))
else:
    g1="NA"
    msq=utils.getMostSevereConsequence(msqs)
    for g in g2c:
        if msq in g2c[g]:
            g1=g
            break
    print("%s\t%s\t%s" %(rsID,g1,msq))
