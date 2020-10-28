#!/usr/bin/python3

import sys
import pandas as pd
import os

if len(sys.argv[1:])==0:
    print("No arguments provided",file=sys.stderr)
    sys.exit(0)

dfs=list()
for i in range(1,len(sys.argv)):
    #print(sys.argv[i],file=sys.stderr)
    bname=os.path.basename(sys.argv[i])
    name=bname[:-4] # name of the experiment
    df=pd.read_table(sys.argv[i],header=0,comment="#",low_memory=False)
    df.set_index(keys="Gene ID",inplace=True,verify_integrity=True)
    df["Experiment"]=name
    dfs.append(df)
    
DF=pd.concat(dfs)
DF.to_csv(path_or_buf=sys.stdout,sep="\t")

