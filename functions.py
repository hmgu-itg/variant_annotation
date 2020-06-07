import requests
import datetime
import sys
import json 
import re
from requests.exceptions import Timeout,TooManyRedirects,RequestException

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

def makeRSQueryURL(rsID,build="38"):
    ext = "/variant_recoder/homo_sapiens/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext+rsID+"?"

def makeRSListQueryURL(build="38"):
    ext = "/variant_recoder/homo_sapiens/"
    server = "http://grch"+build+".rest.ensembl.org"
    if build=="38":
        server = "http://rest.ensembl.org"

    return server+ext

# this function is primarily used to get variant's position
def parseSPDI(string):
    L=string.rsplit(":")
    c=L[0]
    m=re.search("NC_0+(\d+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=int(L[1])
    ref=L[2]
    alt=L[3]
    lref=len(ref)
    lalt=len(alt)

    # ref can be the length of the deleted sequence
    m=re.search("^(\d+)$",ref)
    if m:
        if lalt==1: # SNP
            pos=pos+1
    elif lref==1 and lalt==1:# SNP
        pos=pos+1

    return {"chr":c,"pos":pos,"ref":ref,"alt":alt}

def restQuery(URL,data=None,qtype="get",timeout=None):
    func=None

    if qtype=="get":
        func=requests.get
    elif qtype=="post":
        func=requests.post
    else:
        print(str(datetime.datetime.now())+" : getQuery: query type ("+qtype+") has to be either \"get\" or \"post\"",file=sys.stderr)
        sys.stderr.flush()
        return None

    r=None
    try:
        if qtype=="get":
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},timeout=timeout)
        else:
            if not data:
                print(str(datetime.datetime.now())+" : getQuery: Error: postquery requires data",file=sys.stderr)
                sys.stderr.flush()
                return None                
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},data=data,timeout=timeout)

        if not r.ok:
            print(str(datetime.datetime.now())+" : getQuery: Error "+str(r.status_code)+" occured",file=sys.stderr)
            sys.stderr.flush()
            return None

        try:
            ret=r.json()
            return ret
        except ValueError:
            print(str(datetime.datetime.now())+" : getQuery: JSON decoding error", file=sys.stderr)
            sys.stderr.flush()
            return None

    except Timeout as ex:
        print(str(datetime.datetime.now())+" : getQuery: Timeout exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except TooManyRedirects as ex:
        print(str(datetime.datetime.now())+" : getQuery: TooManyRedirects exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except RequestException as ex:
        print(str(datetime.datetime.now())+" : getQuery: RequestException occured", file=sys.stderr)
        sys.stderr.flush()
        return None

# check if provided input is a valid variant ID
# valid ID: either rsID or chr_pos_ref_alt, where at least one of ref and alt has length 1
def checkID(id):
    m=re.search("^rs\d+",id)
    if m:
        return True
    
    m=re.search("^\d+:\d+_([ATGC]+)_([ATGC]+)",id)
    if m:
        if (len(m.group(1))==1 or len(m.group(2))==1):
            return True
        else:
            return False
    else:
        return False

