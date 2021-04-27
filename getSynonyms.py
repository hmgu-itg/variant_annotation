#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging

from varannot import variant
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser(description="Get synonyms for a given rs ID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="varID",required=True)

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

rs=args.id
logging.getLogger("variant").setLevel(logging.DEBUG)

#---------------------------------------------------------------------------------------------------------------------------

data=query.restQuery(query.makeRsPhenotypeQuery2URL(rs,build))

L=list()
if data:
    if "synonyms" in data:
        L=list(filter(lambda x:x!=rs,data["synonyms"]))

for x in L:
    print(x)
