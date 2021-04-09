#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd

#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Select signal peaks from association results")
parser.add_argument('--input','-i',type=str,action="store",help="Input file with association results",required=True)
parser.add_argument('--threshold','-t',type=float,action="store",help="P-value threshold",required=True)
parser.add_argument('--flank','-f',type=int,action="store",help="flank size (bp)",required=True)
parser.add_argument('--chr','-c',type=str,action="store",help="Chromosome label column",required=True)
parser.add_argument('--pos','-p',type=str,action="store",help="Variant position column",required=True)
parser.add_argument('--pval','-v',type=str,action="store",help="P-value column",required=True)
parser.add_argument('--output','-o',type=str,action="store",help="Output directory; default: same directory as the input file",required=False)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

infile=args.input
pt=args.threshold
flank=args.flank
chrname=args.chr
posname=args.pos
pvalname=args.pval
if args.output is None:
    outdir=os.path.dirname(infile)
else:
    outdir=args.output

#----------------------------------------------------------------------------------------------------------------------------------

def find_peaks(df,flank,posname,pvalname):
    # df is position sorted
    L=[]
    peak_p=1.0
    peak=None

    for idx,row in df.iterrows():
        if peak is None:
            peak=row[posname]
            peak_p=row[pvalname]
        else:
            if peak+flank<row[posname]:
                L.append(peak)
                peak=row[posname]
                peak_p=row[pvalname]
            else:
                if row[pvalname]<peak_p:
                    peak=row[posname]
                    peak_p=row[pvalname]
    L.append(peak)
    criterion=df[posname].map(lambda x: x in L)
    ret=df.loc[criterion].copy()
    return ret

chunks=[]
for c in pd.read_csv(fname,sep="\t",iterator=True,chunksize=10000):
    t=c.loc[pd.to_numeric(c[pvalname],errors="coerce").notnull()]
    t=t.loc[t[pvalname].astype(float)<pt]
    chunks.append(t)
    
df=pd.concat(chunks)
df=df.astype({pvalname:float,posname:int})    

for c,tab in df.groupby(chrname):
    r=find_peaks(tab.sort_values(by=posname),flank,posname,pvalname)
    print(r.to_string())
    r.to_csv(os.path.join(outdir,"peaks_chr"+str(c)+".txt"),sep="\t",index=False)

sys.exit(0)
