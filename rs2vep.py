#!/usr/bin/python3

import sys
import os
import argparse
import json
import logging
import pandas as pd
import re

from varannot import utils
from varannot import query
from varannot import config

#----------------------------------------------------------------------------------------------------------------------------------

def join_consequences(string):
    # print("String: "+string)
    # ",".join(list(map(lambda x:x.replace('\'',''),re.split(',\s*',re.sub(r'[][]','',string)))))
    return ",".join(string)

def main():
    build="38"
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Get VEP consequences for rsID")
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38",required=False)
    parser.add_argument('--id','-i', action="store",help="rs ID",required=False)
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")

    try:
        args=parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if args.build!=None:
        build=args.build

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    rsID=args.id

    LOGGER=logging.getLogger("rs2vep")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.utils").addHandler(ch)
    logging.getLogger("varannot.utils").setLevel(verbosity)

    #---------------------------------------------------------------------------------------------------------------------------

    if sys.stdin.isatty():
        r=query.restQuery(query.makeVepRSQueryURL(rsID,build=build))
        LOGGER.debug("Record length: %d" %(len(r)))
        LOGGER.debug(json.dumps(r,indent=4,sort_keys=True))
        for record in r:
            if "transcript_consequences" not in record:
                continue
            df=pd.json_normalize(record["transcript_consequences"])
            df["rsID"]=rsID
            for c in ["rsID","gene_symbol","gene_id","transcript_id","variant_allele","sift_prediction","sift_score","polyphen_prediction","polyphen_score","consequence_terms"]:
                if c not in df:
                    df[c]="NA"
            df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["rsID","gene_symbol","gene_id","transcript_id","consequence_terms","variant_allele","sift_prediction","sift_score","polyphen_prediction","polyphen_score"])
    else:
        dfs=list()
        i=1
        for L in utils.chunks([line.rstrip() for line in sys.stdin.readlines()],config.VEP_POST_MAX):
            LOGGER.info("Current chunk: %d" % i)
            i+=1
            r=query.restQuery(query.makeVepRSListQueryURL(build=build),data=utils.list2string(L),qtype="post")
            if r is None:
                continue
            LOGGER.debug("Record length: %d" %(len(r)))
            LOGGER.debug("\n======= VEP query ========\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for record in r:
                if "transcript_consequences" not in record:
                    continue
                df=pd.json_normalize(record["transcript_consequences"])
                df["rsID"]=record["id"]
                for c in ["rsID","gene_symbol","gene_id","transcript_id","variant_allele","sift_prediction","sift_score","polyphen_prediction","polyphen_score","consequence_terms"]:
                    if c not in df:
                        df[c]="NA"
                dfs.append(df)
        DF=pd.concat(dfs)
        DF["consequence_terms"]=DF["consequence_terms"].transform(lambda x:",".join(x))
        DF.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["rsID","gene_symbol","gene_id","transcript_id","consequence_terms","variant_allele","sift_prediction","sift_score","polyphen_prediction","polyphen_score"])

if __name__=="__main__":
    main()
