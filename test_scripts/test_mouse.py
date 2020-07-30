#!/usr/bin/python3

import sys
import argparse
import requests
import json

from modules import mouse

#----------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Test mouse models")
requiredArgs=parser.add_argument_group('required arguments')
requiredArgs.add_argument('--id','-i', action="store",help="gene ID",required=True)

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

ID=args.id

#---------------------------------------------------------------------------------------------------------------------------

# r=query.restQuery(config.EXAC_VAR_URL % "1-1158562-AAC-A")
# if r is not None:
#     print(json.dumps(r, indent=4, sort_keys=True))

mouseDF=mouse.getMousePhenotypes(ID)
print(mouseDF)

