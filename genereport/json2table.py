#!/usr/bin/python3

import sys
import json
import pandas as pd
from varannot import config

#---------------------------------------------------------------------------------------------------------------------------

def main():
    data=json.load(sys.stdin)
    df=pd.json_normalize(data["data"])
    df["consequence2"]=df["consequence"].map(config.VEP_CONSEQUENCES)
    df.to_csv(sys.stdout,index=False,sep="\t")

if __name__=="__main__":
    main()
