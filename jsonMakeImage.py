#!/usr/bin/python3

import sys
import argparse
import json

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Create JSON image from STDIN")
    parser.add_argument('--caption','-c', action="store",help="Caption",required=False,default=None)

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    caption=args.caption

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------
    
    input_data=json.load(sys.stdin)
    input_data["type"]="image"
    if caption:
       input_data["caption"]=caption
    print(input_data)
    
if __name__=="__main__":
    main()
