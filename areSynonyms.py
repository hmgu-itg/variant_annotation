#!/usr/bin/python3

import sys
#import os
#import argparse
#import re
import logging

from varannot import variant
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"

# parser = argparse.ArgumentParser(description="Returns 0 if the two provided rsIDs are synonyms, 1 otherwise")
# parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
# requiredArgs=parser.add_argument_group('required arguments')
# requiredArgs.add_argument('--first','-f', action="store",help="rs1",required=True)
# requiredArgs.add_argument('--second','-s', action="store",help="rs2",required=True)

# if len(sys.argv[1:])==0:
#     parser.print_help()
#     sys.exit(0)

# try:
#     args=parser.parse_args()
# except:
#     parser.print_help()
#     sys.exit(0)

# if args.build!=None:
#     build=args.build

# rs1=args.first
# rs2=args.second

rs1=sys.argv[1]
rs2=sys.argv[2]

if rs1==rs2:
    sys.exit(0)

LOGGER=logging.getLogger("areSynonyms")
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.variant").addHandler(ch)
logging.getLogger("varannot.variant").setLevel(logging.DEBUG)
logging.getLogger("varannot.query").addHandler(ch)
logging.getLogger("varannot.query").setLevel(logging.DEBUG)

#---------------------------------------------------------------------------------------------------------------------------

if rs1==rs2:
    sys.exit(0)

data1=query.restQuery(query.makeRsPhenotypeQuery2URL(rs1,build))
data2=query.restQuery(query.makeRsPhenotypeQuery2URL(rs2,build))

L1=list()
L2=list()
if data1:
    if "synonyms" in data1:
        L1=list(filter(lambda x:x!=rs1,data1["synonyms"]))

if rs2 in L1:
    sys.exit(0)
        
if data2:
    if "synonyms" in data2:
        L2=list(filter(lambda x:x!=rs2,data2["synonyms"]))

if rs1 in L2:
    sys.exit(0)

I=[x for x in L1 if x in L2]

if len(I)!=0:
    sys.exit(0)
else:
    sys.exit(1)

