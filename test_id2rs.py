#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from functions import *
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
m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varID)
chrom=m.group(1)
pos=int(m.group(2))
a1=m.group(3)
a2=m.group(4)

#---------------------------------------------------------------------------------------------------------------------------
#r=restQuery(makeOverlapVarQueryURL(chrom,pos,pos,build))
#print(json.dumps(r, indent=4, sort_keys=True))

for x in id2rs(varID,build):
#    print("")
#    print("id",x["id"],sep="\t")
#    print("chrom",x["seq_region_name"],sep="\t")
#    print("start",x["start"],sep="\t")
#    print("end",x["end"],sep="\t")
#    print("alleles",x["alleles"],sep="\t")
    print(x["id"],varID,sep="\t");
