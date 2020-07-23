#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from variant import *
#----------------------------------------------------------------------------------------------------------------------------------

build="38"
window=0

parser = argparse.ArgumentParser(description="Get rs ID for given variant ID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
parser.add_argument('--window','-w', action="store",help="bp window around variant: default: 0", default="0")
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

if args.window:
    window=int(args.window)

varID=args.id
m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)",varID)
chrom=m.group(1)
pos=int(m.group(2))
a1=m.group(3)
a2=m.group(4)

#---------------------------------------------------------------------------------------------------------------------------
#r=restQuery(makeOverlapVarQueryURL(chrom,pos-window,pos+window,build))
#print(json.dumps(r, indent=4, sort_keys=True))

for x in id2rs(varID,build):
    print(x,varID,sep="\t",file=sys.stdout)

# for x in stringdiff("AAAGAAAAGAAA","AAA"):
#     print(x)

# print("")
# for x in stringdiff("AGCCA","AGCCAGCCA"):
#     print(x)


