#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging

from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"
w=0

parser = argparse.ArgumentParser(description="Get reference bases at a given position and a window around")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
parser.add_argument('--window','-w', action="store",help="bp window: default: 0",type=int,default=0)
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--chrom','-c', action="store",help="chrom",required=True)
requiredArgs.add_argument('--pos','-p', action="store",help="pos",type=int,required=True)

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

w=args.window
c=args.chrom
p=args.pos

#---------------------------------------------------------------------------------------------------------------------------
if w!=0:
    s1=query.getRefSeq(c,p-w,p-1,build=build)
    s2=query.getRefSeq(c,p,p,build=build)
    s3=query.getRefSeq(c,p+1,p+w,build=build)
    print(s1,s2,s3,sep=" ")
else:
    print(query.getRefSeq(c,p,p,build=build))
