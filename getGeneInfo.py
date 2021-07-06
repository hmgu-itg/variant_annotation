#!/usr/bin/python3

import sys
import os
import argparse
import json
import logging
import re

from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

build="38"
verbosity=logging.INFO

parser = argparse.ArgumentParser(description="For a  gene ID/name,get basic info")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
parser.add_argument("--input", "-i", help="Input gene ID/gene name",required=True,action='store')

try:
    args=parser.parse_args()
except:
    sys.exit(1)

if args.build!=None:
    build=args.build

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("getHomologs")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler()
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("varannot.query").addHandler(ch)
logging.getLogger("varannot.query").setLevel(verbosity)

input_str=args.input
#---------------------------------------------------------------------------------------------------------------------------
# if re.match("ENSMUSG\d+",input_str):
#     d=query.restQuery(query.makeHomologyURL(input_str,build=build,species="human",homology_type="all"))
#     if d is None:
#         print("%s\t%s\tNA\tNA" %(input_str,input_str))
#         sys.exit(0)       
#     for h in d["data"]:
#         for x in h["homologies"]:
#             print("%s\t%s\t%s\t%s" %(input_str,input_str,x["target"]["id"],x["type"]))
#     print(json.dumps(d,indent=4,sort_keys=True))
# else:
#     d=makeGeneNameXQueryURL(name,build="38",species="homo_sapiens"):
#     d=query.restQuery(query.makeHomologySymbolURL(input_str,source_species="mouse",build=build,target_species="human",homology_type="all"))
#     if d is None:
#         print("%s\t%s\tNA\tNA" %(input_str,g["id"]))
#         sys.exit(0)       
#     for h in d["data"]:
#         for x in h["homologies"]:
#             print("%s\t%s\t%s\t%s" %(input_str,h["id"],x["target"]["id"],x["type"]))


data=query.restQuery(query.makeGeneSymbolQueryURL(input_str,build=build))
if data is None:
    print("%s\tNA\tNA\tNA\tNA" %(input_str))
    sys.exit(0)
#print(json.dumps(data,indent=4,sort_keys=True))
print("%s\t%s\t%s\t%s\t%s" %(input_str,data["id"],data["seq_region_name"],data["start"],data["end"]))
