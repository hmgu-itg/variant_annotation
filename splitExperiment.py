#!/usr/bin/python3

import sys
import pandas as pd
import os

# for experiments including "developmental stage" or age, columns correspond to comma separated combinations of stage/age and organ (example: E-MTAB-4840)
# this script splits one input file into several, one for each stage/age

if len(sys.argv[1:])!=1:
    print("No argument provided",file=sys.stderr)
    sys.exit(1)

fname=sys.argv[1]
df=pd.read_table(fname,header=0,comment="#",low_memory=False)
ccols=list(filter(lambda x: "," in x,df.columns))

print(fname)
# if no commas in header, just exit
if len(ccols)==0:
    sys.exit(0)

ncols=list(filter(lambda x: not "," in x,df.columns))

# prefixes
pxs=set()

# old column name --> new column name (just organ)
D=dict()
for s in df.columns:
    if s in ccols:
        t=s.rsplit(", ",1)
        pxs.add(t[0])
        D[s]=t[1]
    else:
        D[s]=s

# one output for each prefix (stage/age)
for px in pxs:
    fname1=fname[:-4]+"_"+px+".tsv"
    C=list()
    for c in df.columns:
        if c in ncols or c.startswith(px):
            C.append(c)
    df2=df[C].rename(columns=D)
    df2.to_csv(fname1,sep="\t",index=False)
    print(fname1)
