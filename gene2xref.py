#!/usr/bin/python3

import sys
import os
import argparse
import logging
import json
import pandas as pd

from varannot import gene

class DictAction(argparse.Action):
    def __call__(self,parser,namespace,values,option_string=None):
        value_dict = { "debug":logging.DEBUG,"info":logging.INFO,"warning":logging.WARNING,"error":logging.ERROR}
        setattr(namespace,self.dest,value_dict.get(values))

#----------------------------------------------------------------------------------------------------------------------------------

DBNAMES={"entrez":"EntrezGene","reactome":"Reactome_gene"}

def main():
    parser=argparse.ArgumentParser(description="Cross reference information for gene ID")
    parser.add_argument("--verbose", "-v",action=DictAction,help="Optional: verbosity level",required=False,choices=("debug","info","warning","error"),default=logging.INFO)
    parser.add_argument('--id','-i', action="store",help="Gene ID",required=False,default=None)
    parser.add_argument('--source','-s', action="store",help="Source",required=True,choices=["entrez","reactome"])
    parser.add_argument('--build','-b', action="store",help="Genome build: default: 38",required=False,default="38")

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    build=args.build
    ID=args.id
    verbosity=args.verbose
    source=args.source

    LOGGER=logging.getLogger("gene2xref")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler()
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("varannot.gene").addHandler(ch)
    logging.getLogger("varannot.gene").setLevel(verbosity)

    if sys.stdin.isatty() and ID is None:
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    data=gene.getGeneXrefsDB(ID,DBNAMES[source],build=build)
    if data is None:
        sys.exit(1)
    df=pd.json_normalize(data)
    df["ID"]=ID
    df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=["ID","primary_id","description"])

if __name__=="__main__":
    main()
