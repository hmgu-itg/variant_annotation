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

build="38"
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="Get chr/pos/ref/alt info for rs IDs (read from STDIN)")
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

LOGGER=logging.getLogger("rsList2consequence")
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
i=1
for L in utils.chunks([line.rstrip() for line in sys.stdin.readlines()],config.VEP_POST_MAX):
    LOGGER.info("Current chunk: %d" % i)
    i+=1
    r=variant.rsList2position(L,build,alleles=True)
    if r is None:
        continue
    for ID in r:
        for z in r[ID]:
            print("%s\t%s\t%d\t%s\t%s" %(ID,z["chr"],z["pos"],z["ref"],z["alt"]))        
