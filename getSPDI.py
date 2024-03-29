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

parser = argparse.ArgumentParser(description="Get SPDIs for a given rs")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="rsID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

if args.build!=None:
    build=args.build

varID=args.id

LOGGER=logging.getLogger("rs2position")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.variant").addHandler(ch)
logging.getLogger("varannot.variant").setLevel(logging.DEBUG)

#---------------------------------------------------------------------------------------------------------------------------

for x in variant.rs2spdi(varID,build):
    print(x,file=sys.stdout)
