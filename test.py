#!/usr/bin/python3

import sys, time
import os
import argparse
import re
import datetime
from functions import *
import logging

LOGGER=logging.getLogger("test")
LOGGER.setLevel(logging.DEBUG)
#fh=logging.FileHandler('test.log')
#fh.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
#fh.setFormatter(formatter)
ch.setFormatter(formatter)
#LOGGER.addHandler(fh)
LOGGER.addHandler(ch)

#----------------------------------------------------------------------------------------------------------------------------------

build="38"

parser = argparse.ArgumentParser(description="Get chromosome, position and alleles for given rsID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
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
    
#---------------------------------------------------------------------------------------------------------------------------
# r=restQuery(makeRSQueryURL(rsID,build))

# H={}

# if r:
#     if verbose:
#         print("INFO: "+repr(r))

#     if len(r)>1:
#         print("WARNING: More than 1 hash for "+rsID,file=sys.stderr,flush=True)
        
#     x=r[0]
#     r1=x["id"][0]
#     if r1!=rsID:
#         print("WARNING: INPUT ID="+rsID,"RETRIEVED ID="+r1,file=sys.stderr,flush=True)

#     H[rsID]=[]
#     if "spdi" in x:
#         spdi=x["spdi"]
#         for z in spdi:
#             m=re.search("^NC_0+",z)
#             if m:
#                 p=parseSPDI(z)
#                 H[rsID].append(p)

#     s=H[rsID]
#     positions=set(x["chr"]+":"+str(x["pos"]) for x in s)
#     if len(positions)>1:
#         print("ERROR: more than one position for "+rsID,file=sys.stderr,flush=True)
#     elif len(positions)<1:
#         print("ERROR: no position for "+rsID,file=sys.stderr,flush=True)
#     else:
#         L=positions.pop().rsplit(":")
#         print(rsID,L[0],L[1],sep='\t',file=sys.stdout,flush=True)
# else:
#     print("ERROR: restQuery returned None for "+rsID,file=sys.stderr,flush=True)    
#     print(rsID,"NA","NA",sep='\t')

#--------------------------------------------------------------------------------------------------------------

# chrom="1"
# start=int(1000000)
# end=int(1015000)
# URL=makePhenoOverlapQueryURL(chrom,start,end,build=build)
# print(URL)
# variants=restQuery(URL,qtype="get")
# if variants:
#     #print(json.dumps(variants, indent=4, sort_keys=True))
#     if len(variants) == 0: 
#         print(str(datetime.datetime.now())+" : getVariantsWithPhenotypes: No variants with phenotypes were found in the region", file=sys.stderr)

#     rsIDs = []
#     for var in variants:
#         if "id" in var:
#             rsIDs.append(var["id"])
#         else:
#             print(var)

#     r = restQuery(makeRSPhenotypeQueryURL(build=build),data=list2string(rsIDs),qtype="post")
#     print(json.dumps(r, indent=4, sort_keys=True))

#print(getRefSeq(1,1000000,1000100))

#--------------------------------------------------------------------------------------------------------------

# q=makeRsPhenotypeQuery2URL(rsID,build)
# print("Query:",q,sep="\t")
# r=restQuery(q)

# print(json.dumps(r,indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

# z=restQuery(makeRSQueryURL(rsID,build=build))
# for x in z:
#     spdis=x["spdi"]
#     for spdi in spdis:
#         #print("SPDI: "+spdi)
#         h=parseSPDI(spdi,alleles=True)
#         ref=h["ref"]
#         alt=h["alt"]
#         p=h["pos"]
#         c=h["chr"]
#         print("%s\t%d\t%s\t%s\t%s" %(c,p,rsID,ref,alt))

#--------------------------------------------------------------------------------------------------------------

# D=parseGTEx("/home/andrei/out.bed.gz","1",14409,29553,"ENSG00000227232")
# print(D)
# D=parseGTEx("/home/andrei/out.bed.gz","1",108500,108600,"chr1_108506_C_T_b38")
# print(D)

#--------------------------------------------------------------------------------------------------------------

#x=getUniprotData("O15169")
#print(x)

#--------------------------------------------------------------------------------------------------------------

# L=getGwasHits("20",3752792)
# for x in L:
#     print(x)

#--------------------------------------------------------------------------------------------------------------

# L=rs2position(rsID)
# if L:
#     for x in L:
#         print("%s\t%s\t%d" %(rsID,x["chr"],x["pos"]))

#--------------------------------------------------------------------------------------------------------------

#print(json.dumps(getGeneInfo("ENSG00000225972"),indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

#print(json.dumps(getGeneXrefs("ENSG00000185973"),indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

#print(getMouseID("ENSG00000185973"))

#--------------------------------------------------------------------------------------------------------------

#print(getMgiID("ENSMUSG00000079834"))

#--------------------------------------------------------------------------------------------------------------

#print(getMgiPhenotypes("MGI:2180203"))

#--------------------------------------------------------------------------------------------------------------
#LOGGER.info("Calling getApprisInfo")

#print(getApprisInfo("ENSG00000185973"))

#--------------------------------------------------------------------------------------------------------------
#LOGGER.info("Calling getVariantInfo")

#info=getVariantInfo(rsID)
#print(json.dumps(info,indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

#print(json.dumps(getPubmed(rsID),indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------
# LOGGER.info("Calling getGeneList")

# x=getGeneList("1",1230000,1250000)
# print(json.dumps(x,indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

#data=restQuery(makeRsPhenotypeQuery2URL(rsID,build))
#print(json.dumps(data,indent=4,sort_keys=True))

#--------------------------------------------------------------------------------------------------------------

data=getExacAF("1",13528,"C","T")
print(json.dumps(data,indent=4,sort_keys=True))
