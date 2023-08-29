#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd
import re

from varannot import mouse
from varannot import utils

#----------------------------------------------------------------------------------------------------------------------------------

def find_convert_links(string):
    if pd.isna(string):
        return string
    regex=re.compile('(DOID:\d+)')
    return regex.sub(lambda m: utils.makeLink("https://www.informatics.jax.org/disease/"+m.group(),m.group()),string)

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="MGI mouse phenotypes")
    parser.add_argument("--verbose", "-v", help="Optional: verbosity level", required=False,choices=("debug","info","warning","error"),default="info")
    parser.add_argument('--json','-j', action="store_true",help="Output JSON",required=False)
    parser.add_argument('--id','-i', action="store",help="Gene ID",required=False,default=None)
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38",required=False,default="38")
    parser.add_argument('--link','-l', action="store_true",help="Output links",required=False)

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
    out_json=args.json
    link=args.link

    LOGGER=logging.getLogger("gene2mouse")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.mouse").addHandler(ch)
    logging.getLogger("varannot.mouse").setLevel(verbosity)
    logging.getLogger("varannot.utils").addHandler(ch)
    logging.getLogger("varannot.utils").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    data=mouse.getMousePhenotypes(ID,build=build)
    if data is None:
        sys.exit(1)
    data["Phenotypes"]=data["Phenotypes"].replace("\s*\|\s*",", ",regex=True)
    if out_json:
        if link:
            data["Allele ID"]=data["Allele ID"].apply(lambda x:utils.makeLink("https://www.informatics.jax.org/allele/"+x,x))            
            data["mouse gene ID"]=data["mouse gene ID"].apply(lambda x:utils.makeLink("https://www.ensembl.org/Mus_musculus/Gene/Summary?g="+x,x))            
            data["MGI ID"]=data["MGI ID"].apply(lambda x:utils.makeLink("https://www.informatics.jax.org/marker/"+x,x))            
            data["Human disease"]=data["Human disease"].apply(lambda x:find_convert_links(x))
        data.rename(columns={"mouse gene ID":"Mouse gene ID"},inplace=True)
        data[["Allele ID","Phenotypes","Human disease","Mouse gene ID","MGI ID"]].to_json(sys.stdout,orient="records")
    else:
        data.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["Allele ID","Phenotypes","Human disease","mouse gene ID","MGI ID"])

if __name__=="__main__":
    main()
