#!/usr/bin/python3

import sys
import argparse
import json

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Encode JSON table")
    parser.add_argument('--no-header','-n', action="store_true",help="Don't display table header",required=False,default=False)
    parser.add_argument('--column-option','-o', action="store",help="Column option",required=False,default=None)
    parser.add_argument('--caption','-c', action="store",help="Caption",required=False,default=None)
    parser.add_argument('--title','-t', action="store",help="Section title",required=False,default=None)
    parser.add_argument('--long','-l', action="store_true",help="Long table",required=False,default=False)
    parser.add_argument('--landscape','-a', action="store_true",help="Landscape mode",required=False,default=False)
    parser.add_argument('--newpage','-p', action="store_true",help="Trigger new page",required=False,default=False)
    parser.add_argument('--add','-d',action="append",help="Add key:value pair",required=False,default=[])

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    caption=args.caption
    title=args.title
    no_header=args.no_header
    option=args.column_option
    longtable=args.long
    landscape=args.landscape
    newpage=args.newpage
    to_add=args.add

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------
    
    # input_data=sys.stdin.read().strip()
    try:
        input_data=json.load(sys.stdin)
    except Exception as e:
        print("Error: "+str(e),file=sys.stderr)
        sys.exit(1)
    d={"type":"table","data":input_data}
    if caption:
       d["caption"]=caption
    if title:
       d["title"]=title
    if no_header:
        d["no_header"]=True
    if landscape:
        d["landscape"]=True
    if longtable:
        d["longtable"]=True
    if newpage:
        d["newpage"]=True
    if option:
        d["column_option"]=option
    for x in to_add:
        key,val=x.split(":")
        d[key]=val
    print(json.dumps(d))
    
if __name__=="__main__":
    main()
