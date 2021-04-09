#!/usr/bin/python3

import sys
import pandas as pd

def find_peaks(df,flank):
    # df is position sorted
    L=[]
    peak_p=1.0
    peak=None

    for idx,row in df.iterrows():
        if peak is None:
            peak=row["ps"]
            peak_p=row["P-value"]
        else:
            if peak+flank<row["ps"]:
                L.append(peak)
                peak=row["ps"]
                peak_p=row["P-value"]
            else:
                if row["P-value"]<peak_p:
                    peak=row["ps"]
                    peak_p=row["P-value"]
    L.append(peak)
                    
    criterion=df["ps"].map(lambda x: x in L)
    ret=df.loc[criterion].copy()
    return ret

fname=sys.argv[1]

chunks=[]
pt=0.001
#print("Input: %s" %(fname))
for c in pd.read_table(fname,sep="\t",iterator=True,chunksize=100,compression="gzip"):
    #t=c.loc[lambda c: c["P-value"].astype(float)<0.01]
    t=c.loc[pd.to_numeric(c["P-value"],errors="coerce").notnull()]
    t=t.loc[t["P-value"].astype(float)<pt]
    chunks.append(t)
    
df=pd.concat(chunks)
df=df.astype({"P-value":float,"ps":int})    
#print(df.dtypes)    

#print(df.columns)
#if "P-value" in df.columns:
#    print("\"P-value\" is present")

#print(df[["P-value"]].head())
#print(df.loc[range(100,120)])
nrows=df.shape[0]
ncols=df.shape[1]
#print("Dim: %d x %d" % (nrows,ncols))

flank=50000
for c,tab in df.groupby("#chr"):
    #print(tab.sort_values(by="ps"),end="\n\n")
    r=find_peaks(tab.sort_values(by="ps"),flank)
    #r.reset_index(inplace=True)
    #r["label"]=str(c)+"_"+r.index
    #r["label"]="chr"+str(c)+"_peak"+r.index.map(str)
    #r["label"]=str(c)+"_"+r["label"]
    #r.set_index("label",inplace=True)
    #r.drop("index",axis=1,inplace=True)
    print(r.to_string())
    r.to_csv("peaks_chr"+str(c)+".txt",sep="\t",index=False)
    #print(r.dtypes)    
    #print("------------------------------------------------\n\n")
    
#print(df.sort_values(by="P-value"))
#print(df.groupby(by="#chr"))

# for i in range(len(df.sort_values(by="P-value"))):
#     print(df.loc[i,"#chr"], df.loc[i,"ps"])
    
exit(0)
