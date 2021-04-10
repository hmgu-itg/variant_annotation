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
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="Get rs ID for given variant ID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
parser.add_argument('--id','-i', action="store",help="varID",required=True)
parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")

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

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("id2rs")
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

LOGGER.debug("Current variant: %s" % varID)

rsIDs=list(variant.id2rs_mod2(varID,build))
if len(rsIDs)==0:
    print("NA",varID,sep="\t",file=sys.stdout)
elif len(rsIDs)==1:
    print(rsIDs[0],varID,sep="\t",file=sys.stdout)
else:
    LOGGER.warning("Several rsIDs for %s" % varID)
    for x in rsIDs:
        print(x,varID,sep="\t",file=sys.stdout)
