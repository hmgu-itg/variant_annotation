#!/usr/bin/python3

import sys
#import os
#import argparse
import json
#import logging
#import re
import pandas as pd

from varannot import config

#----------------------------------------------------------------------------------------------------------------------------------

#build="38"
#verbosity=logging.INFO

# parser = argparse.ArgumentParser(description="For a gene ID/name, get basic info")
# parser.add_argument("--input", "-i", help="Input gene ID/gene name",required=True,action='store')

# try:
#     args=parser.parse_args()
# except:
#     sys.exit(1)

# if args.build!=None:
#     build=args.build

# if args.verbose is not None:
#     if args.verbose=="debug":
#         verbosity=logging.DEBUG
#     elif args.verbose=="warning":
#         verbosity=logging.WARNING
#     elif args.verbose=="error":
#         verbosity=logging.ERROR

# LOGGER=logging.getLogger("getGeneInfo")
# LOGGER.setLevel(verbosity)
# ch=logging.StreamHandler()
# ch.setLevel(verbosity)
# formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
# ch.setFormatter(formatter)
# LOGGER.addHandler(ch)

# logging.getLogger("varannot.query").addHandler(ch)
# logging.getLogger("varannot.query").setLevel(verbosity)

# input_str=args.input

#---------------------------------------------------------------------------------------------------------------------------

data=json.load(sys.stdin)
df=pd.json_normalize(data["data"])
df["consequence2"]=df["consequence"].map(config.VEP_CONSEQUENCES)
# df.replace({"consequence":config.VEP_CONSEQUENCES},inplace=True)
df2=df[df["consequence2"]==df["consequence2"].max()]
df2[["varId","consequence","pValue","phenotype","dbSNP","maf"]].to_csv(sys.stdout,sep="\t",index=False)
idx=df.groupby("phenotype")["pValue"].idxmin()
print(df.loc[idx][["varId","dbSNP","consequence","pValue","phenotype"]].sort_values(by="pValue").reset_index())
# print(df[["varId","consequence","pValue","phenotype"]])
