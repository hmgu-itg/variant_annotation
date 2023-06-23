#!/usr/bin/python3

import sys
import os
import argparse
import re
import logging
import json

from varannot import variant
from varannot import query

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    build="38"
    verbosity=logging.INFO
    parser = argparse.ArgumentParser(description="Get chr, pos, ref, alt, minor allele, maf for a given rs")
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default="38")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    requiredArgs=parser.add_argument_group('required arguments')
    requiredArgs.add_argument('--id','-i', action="store",help="rsID",required=True)

    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

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

#---------------------------------------------------------------------------------------------------------------------------

    data=query.restQuery(query.makeRsPhenotypeQuery2URL(varID,build,pops=False,phenotypes=False))
    LOGGER.debug("\n%s\n" % json.dumps(data,indent=4,sort_keys=True))
    maf=data["MAF"] if "MAF" in data else "NA"
    mallele=data["minor_allele"] if "minor_allele" in data else "NA"
    z=variant.rs2position(varID,build,alleles=True)
    LOGGER.debug(json.dumps(z,indent=4,sort_keys=True))

    if z:
        for record in z:
            print(varID,record["chr"],record["pos"],record["ref"],record["alt"],mallele,maf,sep="\t",file=sys.stdout)

if __name__=="__main__":
    main()
