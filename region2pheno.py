#!/usr/bin/python3

import sys
import re
import argparse
import logging
import json
import pandas as pd

from varannot import phenotype

PHENO_SOURCE={"cosmic":"COSMIC","clinvar":"ClinVar","ddg2p":"DDG2P","dgva":"DGVa","hgmd":"HGMD-PUBLIC","mim":"MIM morbid","gwas":"NHGRI-EBI GWAS catalog","orphanet":"Orphanet","dbvar":"dbVar"}

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="Get phenotypes associated with given gene")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level",required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--source','-s',action="append",choices=list(PHENO_SOURCE.keys()),help="Optional: source",required=False,default=None)
    parser.add_argument('--region','-r',action="store",help="Region: chr:start-end",required=False,default=None)
    parser.add_argument('--build','-b',action="store",help="Genome build: default: 38",required=False,default="38")

    try:
        args=parser.parse_args()
    except ValueError as e:
        print(str(e))
        parser.print_help()
        sys.exit(1)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    build=args.build
    region=args.region
    sources=args.source
    if sources is None:
        sources=list(PHENO_SOURCE.keys())
    sources=[PHENO_SOURCE[x] for x in sources]

    LOGGER=logging.getLogger("region2pheno")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.phenotype").addHandler(ch)
    logging.getLogger("varannot.phenotype").setLevel(verbosity)

    if sys.stdin.isatty() and region is None:
        parser.print_help()
        sys.exit(1)
        
    m=re.match("(\d+):(\d+)-(\d+)",region)
    if m:
        chrom=m.group(1)
        start=m.group(2)
        end=m.group(3)
    else:
        LOGGER.error("Malformed region")
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    LOGGER.debug("Sources: %s" % (",".join(sources)))
    r=phenotype.byRegion(chrom,int(start),int(end),build=build)
    if r is None:
        sys.exit(1)
    res=list()
    for x in r:
        ID=x["id"] if "id" in x else "NA"
        if "phenotype_associations" in x:
            for p in x["phenotype_associations"]:
                if p["source"] in sources:
                    z=p
                    z["ID"]=ID
                    res.append(z)
    df=pd.json_normalize(res)
    df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA")
            
if __name__=="__main__":
    main()
