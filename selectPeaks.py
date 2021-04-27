#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd

#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Select signal peaks from association results")
parser.add_argument('--input','-i',type=str,action="store",help="Input file with association results",required=True)
parser.add_argument('--threshold','-t',type=float,action="store",help="P-value threshold",required=True)
parser.add_argument('--flank','-f',type=int,action="store",help="Flank size (bp)",required=True)
parser.add_argument('--chr','-c',type=str,action="store",help="Chromosome column name",required=True)
parser.add_argument('--pos','-p',type=str,action="store",help="Variant position column name",required=True)
parser.add_argument('--pval','-v',type=str,action="store",help="P-value column name",required=True)
parser.add_argument('--output','-o',type=str,action="store",help="Output directory; default: same directory as the input file",required=False)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

infile=os.path.realpath(args.input)
pt=args.threshold
flank=args.flank
chrname=args.chr
posname=args.pos
pvalname=args.pval
if args.output is None:
    outdir=os.path.dirname(infile)
else:
    outdir=args.output

print("\n")
print("Input file:\t%s" % infile)
print("Output dir:\t%s" % outdir)
print("Flank:\t%d" % flank)
print("Threshold:\t%.10f" % pt)
print("Chr column:\t%s" % chrname)
print("Pos column:\t%s" % posname)
print("Pvalue column:\t%s\n\n" % pvalname)

if not os.path.exists(infile):
    print("ERROR: input file %s does not exist" % infile)
    sys.exit(1)

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

available_columns=pd.read_csv(infile,nrows=1,sep='\t').columns.tolist()
if not chrname in available_columns:
    print("ERROR: %s is not a valid column name" % chrname)
    print("ERROR: available colmns are: %s" % ",".join(available_columns))
    sys.exit(1)
    
if not posname in available_columns:
    print("ERROR: %s is not a valid column name" % posname)
    print("ERROR: available colmns are: %s" % ",".join(available_columns))
    sys.exit(1)
    
if not pvalname in available_columns:
    print("ERROR: %s is not a valid column name" % pvalname)
    print("ERROR: available colmns are: %s" % ",".join(available_columns))
    sys.exit(1)
    
chunks=[]
# reading in and keeping only records with valid P under threshold
for c in pd.read_csv(infile,sep="\t",iterator=True,chunksize=10000):
    t=c.loc[pd.to_numeric(c[pvalname],errors="coerce").notnull()]
    t=t.loc[t[pvalname].astype(float)<pt]
    chunks.append(t)
    
df=pd.concat(chunks)
df=df.astype({pvalname:float,posname:int})    

if df.shape[0]==0:
    print("No variants after thresholding")
    sys.exit(0)
else:
    print("After P-value filtering there are %d signals\n" % df.shape[0])

for c,tab in df.groupby(chrname):
    r=find_peaks(tab.sort_values(by=posname),flank,posname,pvalname)
    print(r.to_string())
    r.to_csv(os.path.join(outdir,"peaks_chr"+str(c)+".txt"),sep="\t",index=False)

sys.exit(0)
