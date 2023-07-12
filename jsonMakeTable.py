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
    print(json.dumps(d))
    
if __name__=="__main__":
    main()
