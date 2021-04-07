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

parser = argparse.ArgumentParser(description="Get VEP consequences for variant IDs (read from STDIN)")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
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

LOGGER=logging.getLogger("id2consequences")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(verbosity)

def f(ID):
    L=ID.split("_")
    L.insert(2,".")
    return " ".join(L)+" . . ."

#---------------------------------------------------------------------------------------------------------------------------

for L in utils.chunks(sys.stdin.readlines(),config.VEP_POST_MAX):
    r=query.restQuery(query.makeVepListQueryURL(build=build),data="{\"variants\":["+",".join(list(map(lambda x:f(x),L)))+"]}",qtype="post")
    if r:
        print(json.dumps(r,indent=4,sort_keys=True))

