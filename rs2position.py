#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging
import json
import pandas as pd

from varannot import variant
from varannot import query
from varannot import utils
from varannot import config

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    build="38"
    verbosity=logging.INFO
    parser = argparse.ArgumentParser(description="Get chr, pos, ref, alt, minor allele, maf for a given rs")
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--id','-i', action="store",help="rsID",required=False,default=None)

    try:
        args=parser.parse_args()
    except:
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

    varID=args.id

    LOGGER=logging.getLogger("rs2position")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.variant").addHandler(ch)
    logging.getLogger("varannot.variant").setLevel(verbosity)
    logging.getLogger("varannot.query").addHandler(ch)
    logging.getLogger("varannot.query").setLevel(verbosity)

    if sys.stdin.isatty() and varID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    if sys.stdin.isatty():
        maf="NA"
        mallele="NA"
        data=query.restQuery(query.makeRsPhenotypeQuery2URL(varID,build,pops=False,phenotypes=False))
        if data:
            LOGGER.debug("\n%s\n" % json.dumps(data,indent=4,sort_keys=True))
            maf=data["MAF"] if "MAF" in data else "NA"
            mallele=data["minor_allele"] if "minor_allele" in data else "NA"
        z=variant.rs2position(varID,build,alleles=True)

        if z:
            LOGGER.debug(json.dumps(z,indent=4,sort_keys=True))
            DF=pd.DataFrame(columns=["chr","pos","ref","alt","minor_allele","maf"])
            for c in ["chr","pos","ref","alt"]:
                DF[c]=[x[c] for x in z]
            DF["minor_allele"]=mallele
            DF["maf"]=maf
            DF["rsID"]=varID
            DF.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["rsID","chr","pos","ref","alt","minor_allele","maf"])
    else:
        dfs=list()
        i=1
        for L in utils.chunks([line.rstrip() for line in sys.stdin.readlines()],config.VEP_POST_MAX):
            maf=dict()
            mallele=dict()
            LOGGER.info("Current chunk: %d" % i)
            i+=1
            r=query.restQuery(query.makeRSPhenotypeQueryURL(build=build,phenotypes=False),data=utils.list2string(L),qtype="post")
            if r is None:
                continue
            LOGGER.debug("Record length: %d" %(len(r)))
            LOGGER.debug("\n======= query ========\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for rs in r:
                maf[rs]=r[rs]["MAF"] if "MAF" in r[rs] else "NA"
                mallele[rs]=r[rs]["minor_allele"] if "minor_allele" in r[rs] else "NA"
            r=variant.rsList2position(L,build=build,alleles=True)
            if r is None:
                continue
            LOGGER.debug("Record length: %d" %(len(r)))
            LOGGER.debug("\n======= query ========\n%s\n==========================\n" % json.dumps(r,indent=4,sort_keys=True))
            for rs in r:
                df=pd.json_normalize(r[rs])
                if rs in maf:
                    df["maf"]=maf[rs]
                else:
                    df["maf"]="NA"
                if rs in mallele:
                    df["minor_allele"]=mallele[rs]
                else:
                    df["minor_allele"]="NA"
                df["rsID"]=rs
                dfs.append(df)
        DF=pd.concat(dfs)
        DF.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["rsID","chr","pos","ref","alt","minor_allele","maf"])

if __name__=="__main__":
    main()
