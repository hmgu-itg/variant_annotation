#!/usr/bin/python3

import sys
import json
import pandas as pd
from varannot import config

#---------------------------------------------------------------------------------------------------------------------------

def main():
    data=json.load(sys.stdin)
    token=data["continuation"]
    df=pd.json_normalize(data["data"])
    if "consequence" in df.columns:
        df["consequence2"]=df["consequence"].map(config.VEP_CONSEQUENCES)
    df.to_csv(sys.stdout,index=False,sep="\t",na_rep="NA",columns=df.columns.sort_values())
    print(token,file=sys.stderr)

if __name__=="__main__":
    main()
