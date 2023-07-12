#!/usr/bin/python3

import sys
import argparse
import json
import ast

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Create JSON text from STDIN")
    parser.add_argument('--title','-t', action="store",help="Section title",required=False,default=None)
    parser.add_argument('--link','-l', action="store",help="Link",required=False,default=None)
    parser.add_argument('--output','-o', action="store",help="Create text or list",choices=["text","list"],required=False,default="text")

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
    if args.output=="list":
        input_data=ast.literal_eval(input_data)
    d={"type":args.output,"data":input_data}
    if title:
       d["title"]=title
    if link:
       d["link"]=link
    print(json.dumps(d))
    
if __name__=="__main__":
    main()
