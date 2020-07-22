#!/usr/bin/python3

import sys
import argparse
import requests
import json
import query
import config
#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Test ExAC API")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="varID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

varID=args.id

#---------------------------------------------------------------------------------------------------------------------------
r=query.restQuery(config.EXAC_VAR_URL % "14-21853913-T-C")
if r is not None:
    print(json.dumps(r, indent=4, sort_keys=True))
