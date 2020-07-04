#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from functions import *
#----------------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser(description="Testing VEP consequences for given rsID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--rs','-r', action="store",help="rsID",required=True)

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

rsID=args.rs

info=getVariantInfo(rsID)
#print(json.dumps(info,indent=4,sort_keys=True))
c=info["consequence"]
print(c)
for m in info["mappings"]:
    print("%s/%s" %(m["ref"],m["alt"]))
    vep=getVepData(m)
    print(vep[0]["most_severe_consequence"])
print("")
