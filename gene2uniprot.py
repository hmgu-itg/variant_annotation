#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd

from varannot import uniprot
from varannot import gene

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    verbosity=logging.INFO
    parser = argparse.ArgumentParser(description="Uniprot information for a gene")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--id','-i', action="store",help="Gene ID",required=False,default=None)
    parser.add_argument('--add','-a', action="store_true",help="Add gene ID to the output",required=False)
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38",required=False,default="38")
    parser.add_argument('--json','-j', action="store_true",help="Output JSON",required=False)

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    build=args.build
    ID=args.id
    add_id=args.add
    out_json=args.json

    LOGGER=logging.getLogger("gene2uniprot")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.uniprot").addHandler(ch)
    logging.getLogger("varannot.uniprot").setLevel(verbosity)
    logging.getLogger("varannot.gene").addHandler(ch)
    logging.getLogger("varannot.gene").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    xrefs=gene.getGeneXrefs(ID,build=build)
    LOGGER.debug("\n%s\n" % json.dumps(xrefs,indent=4,sort_keys=True))
    uniprotID=xrefs["UniProtKB/Swiss-Prot"][0][0]
    df=uniprot.uniprot2df(uniprot.getUniprotData(uniprotID))
    out_columns=["Field","Value"]
    if add_id:
        df["ID"]=ID
        out_columns=["ID","Field","Value"]    
    if out_json:
        df.fillna(value="NA",inplace=True)
        df[out_columns].to_json(sys.stdout,orient="records")
    else:
        df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=out_columns)

if __name__=="__main__":
    main()
