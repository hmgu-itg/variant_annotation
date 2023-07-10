#!/usr/bin/python3

import sys
import os
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
    parser.add_argument('--id','-i',action="store",help="Gene ID",required=False,default=None)
    parser.add_argument('--build','-b',action="store",help="Genome build: default: 38",required=False,default="38")
    parser.add_argument('--json','-j', action="store_true",help="Output JSON",required=False)

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
    out_json=args.json
    ID=args.id
    sources=args.source
    if sources is None:
        sources=list(PHENO_SOURCE.keys())
    sources=[PHENO_SOURCE[x] for x in sources]

    LOGGER=logging.getLogger("gene2pheno")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.phenotype").addHandler(ch)
    logging.getLogger("varannot.phenotype").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    LOGGER.debug("Sources: %s" % (",".join(sources)))
    r=phenotype.byGene(ID,build=build)
    if r is None:
        sys.exit(1)
    res=list()
    for x in r:
        if x["source"] in sources:
            res.append(x)
    df=pd.json_normalize(res)
    if out_json:
        df.fillna(value="NA",inplace=True)
        df.to_json(sys.stdout,orient="records")
    else:
        df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA")
            
if __name__=="__main__":
    main()
