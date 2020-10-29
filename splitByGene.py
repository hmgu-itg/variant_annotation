#!/usr/bin/python3

import sys
import pandas as pd
import os

if len(sys.argv[1:])!=1:
    sys.exit(1)

fname=sys.argv[1]
dname=os.path.dirname(fname)
df=pd.read_table(fname,header=0,low_memory=False,compression="gzip")
count=1
for gene in df["Gene ID"].unique():
    df2=df.loc[df["Gene ID"]==ID].drop(["Gene ID","Gene Name"],axis=1)
    isn=df2.drop("Experiment",axis=1).isnull()
    df3=df2[~isn.all(axis=1)]
    df3.to_csv(dname+"/"+gene+".tsv.gz",sep="\t",index=False,compression="gzip")
    if count % 1000 == 0:
        print("Done\t%d" % count)
    count+=1
