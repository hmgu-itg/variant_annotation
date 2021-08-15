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

LOGGER=logging.getLogger("rs2gtex")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.utils").addHandler(ch)
logging.getLogger("varannot.utils").setLevel(verbosity)

if not utils.isRS(rsID):
    LOGGER.error("%s is not an rs ID" % rsID)
    sys.exit(1)

#---------------------------------------------------------------------------------------------------------------------------

r=query.restQuery(query.makeGTEXQueryURL(rsID,build=build))
LOGGER.debug(json.dumps(r,indent=4,sort_keys=True))

if r:
    for x in r:
        if x["statistic"]=="p-value":
            print("%s\t%s\t%s\t%s" % (rsID,x["tissue"],x["gene"],x["value"]))
