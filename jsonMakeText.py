#!/usr/bin/python3

import sys
import argparse
import json

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Create JSON text from STDIN")
    parser.add_argument('--title','-t', action="store",help="Section title",required=False,default=None)
    parser.add_argument('--link','-l', action="store",help="Link",required=False,default=None)

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    link=args.link
    title=args.title

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------
    
    input_data=sys.stdin.read().strip()
    d={"type":"text","data":input_data}
    if title:
       d["title"]=title
    if link:
       d["link"]=link
    print(json.dumps(d))
    
if __name__=="__main__":
    main()
