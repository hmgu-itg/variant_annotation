#!/usr/bin/python3

import sys
import os
import argparse
import json
import logging

from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="For a  gene ID/name,get a list of mouse homologs")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
parser.add_argument("--input", "-i", help="Input gene ID/gene name",required=True,action='store')

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

LOGGER=logging.getLogger("getHomologs")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(verbosity)

#---------------------------------------------------------------------------------------------------------------------------
R=dict()
if not args.rs:
    R=variant.id2rs_list([line.rstrip() for line in sys.stdin.readlines()],build=build,skip_non_rs=True,keep_all=False)
else:
    for line in sys.stdin.readlines():
        L=line.rstrip().split("\t")
        R[L[0]]={L[1]}

rs_not_found=[varid for varid in R if R[varid]=={"NA"}]
LOGGER.debug("%d variants with no rs IDs" % len(rs_not_found))
# for variants without rs IDs no phenotype annotations are added
R0=variant.addConsequencesToIDList(rs_not_found,build=build,most_severe_only=True,gene_key="gene_symbol")
R1=variant.addConsequencesToRSList([x for s in list([list(x) for x in R.values()]) for x in s],build=build,most_severe_only=True,gene_key="gene_symbol")
R2=variant.addPhenotypesToRSList([x for s in list([list(x) for x in R.values()]) for x in s],build=build)

# for some input IDs, rs ID might be "NA", in which case output the original ID instead
for v in R:
    if v in rs_not_found:
        (c,p,a1,a2)=v.split("_")
        print("%s\t%s\t%s\t%s" % (v,c+":"+p,json.dumps(R0[v]),"[]"))
    else:
        print("%s\t%s\t%s\t%s" % (v,list(R[v])[0],json.dumps(R1[list(R[v])[0]]),json.dumps(list(R2[list(R[v])[0]]))))
