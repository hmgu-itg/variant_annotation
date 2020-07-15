#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from functions import *
import logging
import config
import gxa
import utils

LOGGER=logging.getLogger("test_gxa")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("gxa").setLevel(logging.DEBUG)
logging.getLogger("utils").setLevel(logging.DEBUG)

#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Test GXA module")
parser.add_argument('--output','-o', action="store",help="Output directory")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="gene ID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

ID=args.id

outdir="./"+ID
if args.output:
    outdir=args.output
    
if outdir.endswith("/"):
    outdir=outdir[:-1]

#---------------------------------------------------------------------------------------------------------------------------

if not utils.createDir(outdir):
    LOGGER.error("Could not create output dir %s" % outdir)

config.OUTPUT_DIR=outdir

gxa.getGxa(ID)
