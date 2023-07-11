#!/usr/bin/python3

import sys
import argparse
import json

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Create JSON table from STDIN")
    parser.add_argument('--caption','-c', action="store",help="Caption",required=False,default=None)
    parser.add_argument('--title','-t', action="store",help="Section title",required=False,default=None)

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    caption=args.caption
    title=args.title

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------
    
    # input_data=sys.stdin.read().strip()
    input_data=json.load(sys.stdin)
    d={"type":"table","data":input_data}
    if caption:
       d["caption"]=caption
    if title:
       d["title"]=title
    print(json.dumps(d))
    
if __name__=="__main__":
    main()
