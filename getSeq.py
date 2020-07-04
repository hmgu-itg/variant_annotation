#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from functions import *
#----------------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser(description="Get genomic interval sequence")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
parser.add_argument('--verbose','-v',default=False,action="store_true",help="verbose output")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--chrom','-c', action="store",help="chromosome",required=True)
requiredArgs.add_argument('--start','-s', action="store",help="start",required=True)
requiredArgs.add_argument('--end','-e', action="store",help="end",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

verbose=args.verbose

if args.build!=None:
    build=args.build

c=args.chrom
s=args.start
e=args.end
    
# ---------------------------------------------------------------------------------------------------------------------------

print(getRefSeq(c,int(s),int(e)))
