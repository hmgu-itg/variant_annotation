#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging
import json

from varannot import variant
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="Get chr:pos:ref:alt for a given rs")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="rsID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

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

varID=args.id

LOGGER=logging.getLogger("rs2position")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.variant").addHandler(ch)
logging.getLogger("varannot.variant").setLevel(verbosity)
logging.getLogger("varannot.query").addHandler(ch)
logging.getLogger("varannot.query").setLevel(verbosity)

#---------------------------------------------------------------------------------------------------------------------------

z=variant.rs2position(varID,build,alleles=True)

if z:
    for x in z:
        print(varID,x["chr"],x["pos"],x["ref"],x["alt"],sep="\t",file=sys.stdout)
