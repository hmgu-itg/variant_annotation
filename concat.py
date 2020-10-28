#!/usr/bin/python3

import sys
import pandas as pd
import os

if len(sys.argv[1:])==0:
    sys.exit(0)

for i in range(1,len(sys.argv)):
    print(sys.argv[i])
    bname=os.path.basename(sys.argv[i])
    name=bname[:-4]
    DF=pd.read_table(sys.argv[i],header=None,comment="#")


