#!/usr/bin/python3

import sys
import argparse
import json

#----------------------------------------------------------------------------------------------------------------------------------

def main():
    parser=argparse.ArgumentParser(description="Create JSON image from STDIN")
    parser.add_argument('--caption','-c', action="store",help="Caption",required=False,default=None)
    parser.add_argument('--title','-t', action="store",help="Section title",required=False,default=None)
    parser.add_argument('--add','-d',action="append",help="Add key:value pair",required=False,default=[])

    try:
        args=parser.parse_args()
    except:
        sys.exit(0)

    caption=args.caption
    title=args.title
    to_add=args.add

    if sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)
        
    #---------------------------------------------------------------------------------------------------------------------------

    try:
        input_data=json.load(sys.stdin)
    except Exception as e:
        print("Error: "+str(e),file=sys.stderr)
        sys.exit(1)
    input_data["type"]="image"
    if caption:
       input_data["caption"]=caption
    if title:
       input_data["title"]=title
    for x in to_add:
        key,val=x.split(":")
        input_data[key]=val
    print(json.dumps(input_data))
    
if __name__=="__main__":
    main()
