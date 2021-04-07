import requests
from requests.exceptions import Timeout,TooManyRedirects,RequestException
import logging
import re

from varannot import config

LOGGER=logging.getLogger(__name__)

def getServerName(build="38"):
    server="https://grch"+build+".rest.ensembl.org"
    if build=="38":
        server="https://rest.ensembl.org"
    return server

def makeOverlapVarQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/"
    return getServerName(build)+ext+"%s:%s-%s?feature=variation"  % (chrom,str(start),str(end))

def makeRefSeqQueryURL(chrom,start,end,build="38"):
    ext="/sequence/region/human/"
    return getServerName(build)+ext+"%s:%s..%s:1?"  % (chrom,str(start),str(end))

def makeRSQueryURL(rsID,build="38"):
    ext="/variant_recoder/homo_sapiens/"
    return getServerName(build)+ext+rsID+"?"

def makeHomologyURL(ID,species="mouse",build="38"):
    ext="/homology/id/%s?&target_species=%s&aligned=0&sequence=none&type=orthologues" %(ID,species)
    return getServerName(build)+ext

def makeGeneQueryURL(ID,build="38"):
    ext="/lookup/id/"
    return getServerName(build)+ext+ID

def makeGeneXQueryURL(ID,build="38"):
    ext="/xrefs/id/"
    return getServerName(build)+ext+"%s?all_levels=1" % ID

def makeGeneXQueryURL2(ID,build="38"):
    ext="/xrefs/id/"
    return getServerName(build)+ext+"%s?external_db=MGI" % ID

def makeRSPhenotypeQueryURL(build="38"):
    ext="/variation/homo_sapiens/"
    return getServerName(build)+ext+"?phenotypes=1"

def makePhenoOverlapQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/"
    return getServerName(build)+ext+"%s:%d-%d?feature=variation;variant_set=ph_variants" %(chrom,start,end)

def makeRsPhenotypeQuery2URL(rs,build="38"):
    ext="/variation/human/"
    return getServerName(build)+ext+"%s?pops=1;phenotypes=1" %rs

def makeRSListQueryURL(build="38"):
    ext="/variant_recoder/homo_sapiens/"
    return getServerName(build)+ext

def makeVepQueryURL(chrom,start,end,allele,strand="1",build="38"):
    ext="/vep/homo_sapiens/region/%s:%s-%s:1/%s?" % (chrom,str(start),str(end),allele)
    return getServerName(build)+ext

def makeVepListQueryURL(build="38"):
    ext="/vep/homo_sapiens/region"
    return getServerName(build)+ext

def makeGeneOverlapQueryURL(chrom,start,end,build="38"):
    ext="/overlap/region/human/%s:%d-%d?feature=gene" %(chrom,start,end)
    return getServerName(build)+ext

# ===========================================================================================================================

# if alleles==True, also return alleles
def parseSPDI(string,alleles=False,build="38"):
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
    isSNP=False
    base=None
    seq=None
    isNum=False

    # ref can be the length of the deleted sequence
    # or empty in case of simple insertions
    m=re.match("^(\d+)$",ref)
    if m:
        isNum=True
        if ref=="1" and lalt==1: # one base deleted, one inserted, i.e. variant is a SNP
            isSNP=True
    elif lref==1 and lalt==1:# SNP
        isSNP=True

    if alleles:
        # deleted sequence
        if isNum:
            delseq=getRefSeq(c,pos+1,pos+int(ref),build)
        else:
            delseq=ref

        if isSNP:
            return {"chr":c,"pos":pos+1,"ref":delseq[0],"alt":alt}
        else:
            # base before the deleted sequence
            base=getRefSeq(c,pos,pos,build)
            # deletion
            if len(delseq)>len(alt):
                if delseq.endswith(alt):
                    delseq2=delseq[0:len(delseq)-len(alt)]
                    return {"chr":c,"pos":pos,"ref":base+delseq2,"alt":base}
                else:
                    LOGGER.warning("parseSPDI: deletion: (%s %d %s %s) INS is not suffix of DEL" %(c,pos,ref,alt))
                    return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
            # insertion
            elif len(delseq)<len(alt):
                if lref==0:
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt}
                elif alt.endswith(delseq):
                    alt2=alt[0:len(alt)-len(ref)]
                    return {"chr":c,"pos":pos,"ref":base,"alt":base+alt2}
                else:
                    LOGGER.warning("parseSPDI: insertion: (%s %d %s %s) DEL is not suffix of INS" %(c,pos,ref,alt))
                    return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
            # indel
            else:
                return {"chr":c,"pos":pos,"ref":base+delseq,"alt":base+alt}
    else:
        if isSNP:
            return {"chr":c,"pos":pos+1,"ref":None,"alt":None}
        else:
            return {"chr":c,"pos":pos,"ref":None,"alt":None}

# ===========================================================================================================================

def restQuery(URL,data=None,qtype="get",timeout=None):
    func=None

    if qtype=="get":
        func=requests.get
    elif qtype=="post":
        func=requests.post
    else:
        LOGGER.error("Query type %s has to be either \"get\" or \"post\"" %(qtype))
        return None

    r=None
    try:
        if qtype=="get":
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},timeout=timeout)
        else:
            if not data:
                LOGGER.error("POST query requires data")
                return None                
            r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},data=data,timeout=timeout)

        if not r.ok:
            LOGGER.error("Error %s occured (input URL: %s)" %(str(r.status_code),URL))
            return None

        try:
            ret=r.json()
            return ret
        except ValueError:
            LOGGER.error("JSON decoding error")
            return None

    except Timeout as ex:
        LOGGER.error("Timeout exception occured")
        return None
    except TooManyRedirects as ex:
        LOGGER.error("TooManyRedirects exception occured")
        return None
    except RequestException as ex:
        LOGGER.error("RequestException occured")
        return None

# ===========================================================================================================================

# return a part of the reference genome
def getRefSeq(chrom,start,end,build="38"):
    r=restQuery(makeRefSeqQueryURL(chrom,start,end,build))
    if r:
        return r["seq"]
    else:
        return None

# ===========================================================================================================================

def fetchGnomAD(jsondata,url=config.GNOMAD_URL):
    response=requests.post(url,json=jsondata,headers={"Content-Type": "application/json"})
    #response=restQuery(url,data=jsondata,qtype="post")
    #if not response:
    #    return None
    json=response.json()
    if "errors" in json:
        LOGGER.error(str(json["errors"]))
        return None
    return json
