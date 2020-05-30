#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
import functions
#----------------------------------------------------------------------------------------------------------------------------------

build="grch38"

parser = argparse.ArgumentParser(description="Get chromosome, position and alleles for given rsID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: grch38", default="grch38")
parser.add_argument('--verbose','-v',default=False,action="store_true",help="verbose output")
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

verbose=args.verbose

if args.build!=None:
    build=args.build

rsID=args.rs
    
ext = "/variant_recoder/homo_sapiens/"
server = "http://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "http://rest.ensembl.org"


#print("Current build: "+build,file=sys.stderr)

#---------------------------------------------------------------------------------------------------------------------------

timeout=60
max_attempts=5

URL=server+ext+rsID+"?"
r=functions.restQuery(URL)
H={}

if r:
    if verbose:
        print("INFO: "+repr(r))

    if len(r)>1:
        print("WARNING: More than 1 hash for "+rsID,file=sys.stderr,flush=True)
        
    x=r[0]
    r1=x["id"][0]
    if r1!=rsID:
        print("WARNING: INPUT ID="+rsID,"RETRIEVED ID="+r1,file=sys.stderr,flush=True)

    H[rsID]=[]
    if "spdi" in x:
        spdi=x["spdi"]
        for z in spdi:
            m=re.search("^NC_0+",z)
            if m:
                p=functions.parseSPDI(z)
                H[rsID].append(p)

    s=H[rsID]
    positions=set(x["chr"]+":"+str(x["pos"]) for x in s)
    if len(positions)>1:
        print("ERROR: more than one position for "+rsID,file=sys.stderr,flush=True)
    elif len(positions)<1:
        print("ERROR: no position for "+rsID,file=sys.stderr,flush=True)
    else:
        L=positions.pop().rsplit(":")
        print(rsID,L[0],L[1],sep='\t',file=sys.stdout,flush=True)
else:
    print("ERROR: getResponse2 returned None for "+rsID,file=sys.stderr,flush=True)    
    print(rsID,"NA","NA",sep='\t')
