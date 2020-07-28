#!/usr/bin/python3

import sys
import os
import argparse
import re

from modules import variant

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
# m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varID)
# chrom=m.group(1)
# pos=int(m.group(2))
# a1=m.group(3)
# a2=m.group(4)

#---------------------------------------------------------------------------------------------------------------------------

for x in variant.id2rs(varID,build):
    print(x,varID,sep="\t",file=sys.stdout)

