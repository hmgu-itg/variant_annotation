import logging
import re
import jinja2

LOGGER=logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)
ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

# check if provided input is a valid variant ID
# valid ID: either rsID or chr_pos_ref_alt
def checkID(id):
    m=re.search("^rs\d+",id)
    if m:
        return True
    
    m=re.search("^\d+_\d+_([ATGC]+)_([ATGC]+)",id)
    if m:
        if len(m.group(1))==1 and len(m.group(2))==1:
            return True

        if m.group(1).startswith(m.group(2)) or m.group(2).startswith(m.group(1)):
            return True
        else:
            return False
    else:
        return False

def chunks(L,n):
    """Yield successive n-sized chunks from L"""
    for i in range(0, len(L), n):
        yield L[i:i + n]

# VEP-style string
def getLocationString(chrom,pos,ref,alt):
    if len(ref)>1 and len(alt)==1:
        start=pos+1
        end=pos+len(ref)-1
        allele1=""
        allele2="-"
    elif len(alt)>1 and len(ref)==1:
        start=pos+1
        end=pos
        allele1=""
        allele2=alt[1:]
    elif len(ref)==1 and len(alt)==1:
        start=pos
        end=pos
        allele1=ref
        allele2=alt
    else:
        LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt),file=sys.stderr)
        return None

    return chrom+":"+str(start)+"-"+str(end)+":"+allele1+"/"+allele2

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

# ----------------------------------------------------------------------------------------------------------------------

def generateHTML(templateFile, data):
    templateLoader = jinja2.FileSystemLoader(searchpath="/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(templateFile)
    return template.render(data)

