#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging

from modules import variant
from modules import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser(description="Get rs ID for given variant ID")
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

varID=args.id

LOGGER=logging.getLogger("id2rs")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("variant").setLevel(logging.DEBUG)

#---------------------------------------------------------------------------------------------------------------------------

rsIDs=list(variant.id2rs(varID,build))
if len(rsIDs)==0:
    print(varID,"NA",sep="\t",file=sys.stdout)
elif len(rsIDs)==1:
    print(varID,rsIDs[0],sep="\t",file=sys.stdout)
else:
    LOGGER.warning("Several rsIDs for %s" % varID)
    for x in rsIDs:
        print(varID,x,sep="\t",file=sys.stdout)
