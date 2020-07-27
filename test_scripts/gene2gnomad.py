#!/usr/bin/env python3

import argparse
import requests
import json
import sys

def get_variant_list(chrom,start,end, dataset="gnomad_r3"):
    fmt_graphql = """
    {
    region(start:%d,stop:%d,chrom:"%s",reference_genome:GRCh38) {
    variants (dataset:%s) {
    variantId
    rsid  
    }
    } 
    }
    """
    req_variantlist = {
        "query": fmt_graphql % (start,end,chrom, dataset)
        }
    response = fetchGnomAD(req_variantlist)
    return response

def fetchGnomAD(jsondata,url="https://gnomad.broadinstitute.org/api"):
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        return None
    return json

# ====================================================================================================

parser=argparse.ArgumentParser(description="For a given region, return all variants")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--chrom','-c', action="store",help="chrom",required=True)
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

chrom=args.chrom
start=int(args.start)
end=int(args.end)

# ====================================================================================================

#print(json.dumps(get_variant_list(chrom,start,end),indent=4,sort_keys=True))
v=get_variant_list(chrom,start,end)["data"]["region"]["variants"]
for x in v:
    print(x["rsid"],x["variantId"],sep="\t")


