import logging
import re
import jinja2
import os
import sys
import pandas as pd
import tempfile as tf
import plotly.graph_objects as go
import subprocess

from shutil import which

from varannot import config
from varannot import query

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def runLiftOver(input_data,build="38"):
    '''
    Run liftOver

    Input: list of dictionaries with keys chr,start,end,id
    Input: build (37 or 38): build to lift over FROM
    Output: lifted over list of dictionaries with keys chr,start,end,id
    '''

    if build!="38" and build!="37":
        LOGGER.error("provided build: %s; build should be either 37 or 38" % build)
        return None
    
    L=list()

    chain="/usr/local/bin/hg38ToHg19.over.chain.gz"
    if build=="37":
        chain="/usr/local/bin/hg19ToHg38.over.chain.gz"

    in_bed=tf.NamedTemporaryFile(delete=False,mode="w",prefix="annotator_")
    LOGGER.debug("Input bed file: %s" % (in_bed.name))
    out_fname=tf.mktemp(prefix="annotator_")
    LOGGER.debug("Output bed file: %s" % (out_fname))
    unmapped_fname=tf.mktemp(prefix="annotator_")
    LOGGER.debug("Unmapped file: %s" % unmapped_fname)

    for x in input_data:
        in_bed.write("%s\t%d\t%d\t%s\n" %(x["chr"],int(x["start"]),int(x["end"]),x["id"]))
    in_bed.close()

    LOGGER.debug("Input: %d record(s)" % len(input_data))
    
    if which("liftOver") is None:
        LOGGER.error("liftOver was not found in PATH")
        return None
        
    LOGGER.debug("Calling: liftOver %s %s %s %s" %(in_bed.name,chain,out_fname,unmapped_fname))
    cmdline="liftOver %s %s %s %s" %(in_bed.name,chain,out_fname,unmapped_fname)
    subprocess.run(cmdline,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    if not os.path.isfile(out_fname):
        LOGGER.error("liftOver failed to create output file %s" % out_fname)
        return None

    if os.path.getsize(out_fname)==0:
        LOGGER.warning("liftOver produced empty output file %s" % out_fname)
        return None

    count=0
    with open(out_fname) as F:
        for line in F:
            (chrom,start,end,ID)=line.rstrip().split("\t")
            L.append({"chr":chrom,"start":start,"end":end,"id":ID})
            count+=1

    LOGGER.debug("Remapped: %d record(s)" % count)
    unmapped=sum(1 for line in open(unmapped_fname) if not line.startswith("#"))
    LOGGER.debug("Unmapped: %d record(s)" % unmapped)
    
    LOGGER.debug("Removing temporary files")
    if os.path.isfile(in_bed.name):
        os.remove(in_bed.name)
    if os.path.isfile(out_fname):
        os.remove(out_fname)
    if os.path.isfile(unmapped_fname):
        os.remove(unmapped_fname)
        
    return L

# ==============================================================================================================================

def checkDEL(R,build="38"):
    '''
    Check if "del" sequence is in genome at "seq":"pos" position

    Input: dict with keys seq,pos,del,ins
    Output: True/False
    '''

    return query.getRefSeq(R["seq"],R["pos"],R["pos"]+len(R["del"])-1,build)==R["del"]

# ==============================================================================================================================

def checkAlleles(ID,build="38"):
    '''
    Given a variant ID (chr_pos_a1_a2), check if a1 or a2 at chr:pos match genome sequence

    Input: variant ID (chr_pos_a1_a2)
    Output: True/False
    '''

    t=splitID(ID)
    if not t:
        return False
    return query.getRefSeq(t["chr"],t["pos"],t["pos"]+len(t["a1"])-1,build)==t["a1"] or query.getRefSeq(t["chr"],t["pos"],t["pos"]+len(t["a2"])-1,build)==t["a2"]

# ==============================================================================================================================

def getMostSevereConsequence(L):
    '''
    Get the most severe VEP consequence from a list of VEP consequences

    Input: list of VEP consequences
    Output: most severe VEP consequence
    '''
    
    c=0
    cons="NA"
    for x in L:
        if not x in config.VEP_CONSEQUENCES:
            LOGGER.warning("input consequence %s is not in the config.VEP_CONSEQUENCES" %(x))
        else:
            if config.VEP_CONSEQUENCES[x]>c:
                c=config.VEP_CONSEQUENCES[x]
                cons=x
    return cons

# ==============================================================================================================================

def getVarType(R):
    '''
    For a given variant dict with keys seq,pos,del,ins, get variant type

    Input: dict with keys seq,pos,del,ins
    Output: SNP,INS,DEL,INDEL
    '''
    
    if len(R["del"])==1 and len(R["ins"])==1:
        return "SNP"

    if len(R["del"])>len(R["ins"]):
        return "DEL"
    elif len(R["del"])<len(R["ins"]):
        return "INS"
    else:
        return "INDEL"

# ==============================================================================================================================

def var2spdi(R):
    '''
    For a given variant dict with keys seq,pos,del,ins (VCF style), get variant's SPDI string

    Input: dict with keys seq,pos,del,ins
    Output: SPDI string
    '''
    
    t=getVarType(R)
    if t=="SNP":
        return R["seq"]+":"+str(R["pos"]-1)+":"+R["del"]+":"+R["ins"]
    elif t=="DEL":
        return R["seq"]+":"+str(R["pos"])+":"+R["del"][1:]+":"
    elif t=="INS":
        return R["seq"]+":"+str(R["pos"])+"::"+R["ins"][1:]
    else:
        return R["seq"]+":"+str(R["pos"]-1)+":"+R["del"]+":"+R["ins"]

# ==============================================================================================================================

def equivalentVariants(r1,r2,build="38"):
    '''
    Given two variant dicts with keys seq,pos,del,ins, test if both variants result in the same altered sequence

    Input: two variant dicts with keys seq,pos,del,ins
    Output: True/False
    '''
    
    if r1["seq"]!=r2["seq"]:
        return False

    if len(r1["ins"])-len(r1["del"])!=len(r2["ins"])-len(r2["del"]):
        return False
    
    c=r1["seq"]
    
    if r1["pos"]<r2["pos"]:
        R1=r1
        R2=r2
    else:
        R1=r2
        R2=r1
        
    x=R2["pos"]-R1["pos"]
    d1=len(R1["del"])
    d2=len(R2["del"])

    LOGGER.debug("R1: %s" % str(R1))
    LOGGER.debug("R2: %s" % str(R2))
    
    if x+d2>=d1:
        LOGGER.debug("x+d2>=d1")
        if x==0:
            ss1=""
        else:
            ss1=query.getRefSeq(c,R1["pos"],R1["pos"]+x-1,build)
        ss2=R2["ins"]
        str1=ss1+ss2
        LOGGER.debug("str1=%s+%s" % (ss1,ss2))
        
        ss1=R1["ins"]
        if R1["pos"]+d1<=R1["pos"]+d2+x-1:
            ss2=query.getRefSeq(c,R1["pos"]+d1,R1["pos"]+d2+x-1,build)
        else:
            ss2=""
        str2=ss1+ss2
        LOGGER.debug("str2=%s+%s" % (ss1,ss2))
        return str1==str2
    else:
        LOGGER.debug("x+d2<d1")
        if x==0:
            ss1=""
        else:
            ss1=query.getRefSeq(c,R1["pos"],R1["pos"]+x-1,build)
        ss2=R2["ins"]
        if R1["pos"]+x+d2<=R1["pos"]+d1-1:
            ss3=query.getRefSeq(c,R1["pos"]+x+d2,R1["pos"]+d1-1,build)
        else:
            ss3=""
        str1=ss1+ss2+ss3
        LOGGER.debug("str1=%s+%s+%s" % (ss1,ss2,ss3))
        
        str2=R1["ins"]
        LOGGER.debug("str2=%s" % str2)
        return str1==str2

# ==============================================================================================================================

def convertSPDI(spdi,build="38"):
    '''
    Convert SPDI string to a dict with keys seq,pos,del,ins

    Input: SPDI string
    Output: dict with keys seq,pos,del,ins (pos is 1-based)
    '''
    
    L=spdi.rsplit(":")
    c=L[0]
    m=re.search("NC_0+(\d+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=L[1]
    D=L[2]
    I=L[3]

    m=re.match("^(\d+)$",D)
    if m:
        return {"seq":c,"pos":str(int(pos)+1),"del":query.getRefSeq(c,int(pos)+1,int(pos)+int(m.group(1)),build),"ins":I}
    else:
        return {"seq":c,"pos":str(int(pos)+1),"del":D,"ins":I}

# ==============================================================================================================================

def convertVariantID(varid,reverse=False):
    '''
    Convert a given variant ID (1_12345_A1_A2) to a dict with keys seq,pos,del,ins

    Input: variant ID (1_12345_A1_A2), reverse: True/False
    Output: dict with keys seq,pos,del,ins; if reverse==True, the A2 is del, A1 is ins
    '''
    
    L=varid.rsplit("_")
    if reverse:
        return {"seq":L[0],"pos":int(L[1]),"del":L[3],"ins":L[2]}
    else:
        return {"seq":L[0],"pos":int(L[1]),"del":L[2],"ins":L[3]}

# ==============================================================================================================================

# input: dict with keys: "seq", "pos", "del", "ins"
# output: string "1 12345 1_12345_AC_A AC A ..."
# WARNING: has not been properly tested

def variant2vep(variant,reverse=False):
    chrom=variant["seq"]
    pos=variant["pos"]
    ref=variant["del"]
    alt=variant["ins"]

    varid="_".join([chrom,str(pos),ref,alt])
    return " ".join([chrom,str(pos),varid,ref,alt,". . ."])

# ==============================================================================================================================

def df2svg(df,var):
    '''
    For dataframe with population data (there must be a column "Population"), create a barplot

    Input: dataframe, variant name
    Output: file name of the created barplot
    '''
    
    if df is None:
        return None
    if len(df)==0:
        return None
    data=list()
    pops=list(df["Population"])
    for i in range(1,len(df.columns)):
        data.append(go.Bar(name=df.columns[i],x=pops,y=df.iloc[:,i]))
    fig=go.Figure(data=data)
    fig.update_layout(barmode='stack')
    #out=tf.NamedTemporaryFile(dir=config.OUTPUT_DIR_FIG,suffix=".svg")
    fd,fname=tf.mkstemp(dir=config.OUTPUT_DIR_FIG,suffix="_"+var+".svg")
    LOGGER.debug("SVG: %s" % fname)
    fig.write_image(fname)
    return fname

# ==============================================================================================================================

def df2barchart(df):
    if len(df)==0:
        return None

    data=list()
    pops=list(df["Population"])
    for i in range(1,len(df.columns)):
        data.append(go.Bar(name=df.columns[i],x=pops,y=df.iloc[:,i]))
    
    fig=go.Figure(data=data)
    fig.update_layout(barmode='stack')
    out=tf.NamedTemporaryFile(delete=False,mode="w")
    fig.write_html(out.name,full_html=False)
    out.close()
    with open (out.name, "r") as f:
        data=f.readlines()
    f.close()
    if os.path.isfile(out.name):
        os.remove(out.name)

    return "".join(data)

# ======================================================================================================================

def makeLink(url,text):
    '''
    Create an HTML link

    Input: url, text
    Output: HTML href link
    '''

    return "<a href='"+url+"'>"+text+"</a>"

# ======================================================================================================================

def checkGnomadID(var):
    '''
    Checks if provided ID conforms to GnomeAD rules (12-1234567-ACT-A)

    Input  : variant ID
    Output : True/False
    '''
    m=re.search("([XYM\d]+)-(\d+)-([ACGT]+)-([ACGT]+)",var)
    if not m:
        LOGGER.debug("%s is not GnomeAD variant ID" % var)
        return False
    else:
        if len(m.group(3))==1 and len(m.group(4))==1:
            return True
        else:
            if m.group(3).startswith(m.group(4)) or m.group(4).startswith(m.group(3)):
                return True
            else:
                LOGGER.debug("%s is not GnomeAD variant ID" % var)
                return False

# ======================================================================================================================

def checkID(varid):
    '''
    Check if provided input is a valid variant ID (either rsID or chr_pos_A1_A2)

    Input: string
    Output: True/False
    '''
    
    if isRS(varid):
        return True
    
    m=re.search("^\d+_\d+_[ATGC]+_[ATGC]+$",varid)
    if m:
        return True
    else:
        return False

# ======================================================================================================================
def splitID(ID):
    '''
    Split variant ID into dict with keys chr,pos,a1,a2

    Input: variant ID (21_12345_AC_A)
    Output: dict with keys chr,pos,a1,a2 or None
    '''
    
    m=re.search("^(\d+)_(\d+)_([ATGC]+)_([ATGC]+)$",ID)
    if m:
        return {"chr":m.group(1),"pos":m.group(2),"a1":m.group(3),"a2":m.group(4)}
    else:
        return None

# ======================================================================================================================

def isRS(varid):
    '''
    Check if provided input is an rs ID

    Input: ID string
    Output: True/False
    '''

    m=re.search("^rs\d+$",varid)
    if m:
        return True
    return False

# ======================================================================================================================

def chunks(L,n):
    '''
    Yields successive n-sized chunks from the input list L
    '''

    for i in range(0, len(L), n):
        yield L[i:i + n]

# ======================================================================================================================

def getLocationString(chrom,pos,ref,alt):
    '''
    Input: chrom, pos, ref, alt
    Output: VEP-style string, used in queries like "https://rest.ensembl.org/documentation/info/vep_region_get"

    SNP: (1,12345,A,C) --> 1:12345-12345:A/C
    INS: (1,12345,A,AC) --> 1:12346-12345:/C
    DEL: (1,12345,ACT,A) --> 1:12346-12347:/-
    '''

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
        LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt))
        return None

    return chrom+":"+str(start)+"-"+str(end)+":"+allele1+"/"+allele2

# ======================================================================================================================

def list2string(varIDs):
    '''
    Convert a list of variant IDs to string

    Input: list of variant IDs
    Output: string {"ids":["rs1","rs2","rs3"]}
    '''
    
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",varIDs)))+"]}"

# ======================================================================================================================

def generateHTML(templateFile,data):
    '''
    For a given HTML template and data, generate HTML page
    '''
    
    templateEnv=jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"))
    template=templateEnv.get_template(templateFile)
    return template.render(data)

# ======================================================================================================================

def checkFiles(files):
    '''
    Check if files exist

    Input: list of file names
    Output: True if all files exist, False otherwise
    '''

    f=True
    for fn in files:
        if not os.path.isfile(fn):
            LOGGER.warning("File %s does not exist" % fn)
            f=False
    return f

# ======================================================================================================================

def createDir(path):
    '''
    Create directory

    Input: directory name
    Output: directory name if successfull or directory already exists, None otherwise
    '''
    
    if os.path.exists(path):
        LOGGER.info("%s already exists" % path)
        return path
    try:
        os.mkdir(path)
    except OSError:
        LOGGER.error("Creation of the directory %s failed" % path)
        return None
    else:
        return path

# ======================================================================================================================

# def getGxaURL(ID):
#     return config.GXA_URL_PREFIX+ID+config.GXA_URL_SUFFIX

# ======================================================================================================================

def generateTemplate(mapping_names,gene_names,fname):
    f = open(fname,"w")
    f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n")
    f.write("<style>\n\nth{ text-align: left; }\n")
    f.write("#common {\nborder-collapse: collapse;\nwidth: 100%;\ncolumn-width=auto;\ncolor: black;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 15px;\n}\n\n")
    f.write("h1 {\nbackground-color: #8b9bc1;\nmargin: 0;\npadding: 4px 8px 2px 24px;\n-webkit-border-radius: 8px 0 8px 0;\nline-height: 1em;\ndisplay: block;\n-webkit-margin-before: 0px;\n-webkit-margin-after: 0px;-webkit-margin-start: 0px;\n-webkit-margin-end: 0px;\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 25px;\nfont-weight: bold;\n}\n\n")
    f.write(".tab {\noverflow: hidden;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\n}\n")
    f.write(".tab button {\nbackground-color: inherit;\nfloat: left;\nborder: none;\noutline: none;\ncursor: pointer;\npadding: 14px 16px;\ntransition: 0.3s;\n}\n")
    f.write(".tab button:hover {\nbackground-color: #ddd;\n}\n")
    f.write(".tab button.active {\nbackground-color: #ccc;\n}\n")
    f.write(".tabcontent {\ndisplay: none;\npadding: 6px 12px;\nborder: 1px solid #ccc;\nborder-top: none;\n}\n")
    f.write(".genetabcontent {\ndisplay: none;\nfloat: left;\npadding: 0px 12px;\nborder: 1px solid #ccc;\nwidth: 98%;\nheight: 98%;\n}\n")
    # collapsible styles
    f.write(".collapsible {\nbackground-color: #33b2ff;\ncolor: white;\ncursor: pointer;\npadding: 18px;\nwidth: 100%;\nborder: none;\ntext-align: center;\noutline: none;\nfont-size: 15px;\nborder-radius: 12px;\n}\n\n")
    f.write(".active, .collapsible:hover {\nbackground-color: #3399ff;\n}\n\n")
    # ------------------
    f.write(".content {\npadding: 0 18px;\ndisplay: none;\nfloat: none;\ntext-align: center;\noverflow: hidden;\nbackground-color: #f1f1f1;\n}\n\n")

    f.write("#space\n{width: 100%;\n height: 5px;\n margin: 0px;\n padding: 0px;\n border: 0px;\nbackground: #FFFFFF;\nclear:both;\n}\n\n") 

    f.write("</style>\n\n")

    f.write("</head>\n<body bgcolor=\"#E6E6FA\" class=body>\n\n")

    f.write("<div class=\"tab\">\n")
    f.write("<button class=\"tablinks2\" onclick=\"openTab2(event, 'Variant')\" id=\"defaultTabOpen\">Variant</button>\n")
    f.write("<button class=\"tablinks2\" onclick=\"openTab2(event, 'Genes')\">Genes</button>\n")
    f.write("</div>\n\n")

    # Variant tab
    f.write("\n<div id=\"Variant\" class=\"tabcontent2\">\n")

    # Mappings
    f.write("<div class=\"tab\">")
    for i in range(0,len(mapping_names)):
        if i==0:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\" id=\"defaultOpen\">%s</button>" % (mapping_names[i],mapping_names[i]))
        else:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\">%s</button>" % (mapping_names[i],mapping_names[i]))
    f.write("</div>\n")

    # Mapping details
    for i in range(0,len(mapping_names)):
        f.write("<div id=\"%s\" class=\"tabcontent\">\n" % mapping_names[i])

        f.write("<button type=\"button\" class=\"collapsible\">Variant Details</button>\n")
        f.write("<div class=\"content\">\n{{ variant_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby variants associated with phenotypes</button>\n")
        f.write("<div class=\"content\">\n{{ phenotype_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby GWAS signals</button>\n")
        f.write("<div class=\"content\">\n{{ gwas_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Predicted consequences</button>\n")
        f.write("<div class=\"content\">\n{{ vep_table%d  }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Population allele frequencies</button>\n")
        f.write("<div class=\"content\">\n{{ population_table%d }}\n</div>\n" %i)
#        f.write("\n{{ population_table%d }}\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">gnomAD allele frequencies</button>\n")
        f.write("<div class=\"content\">\n{{ gnomad_table%d  }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">ENSEMBL Regulation</button>\n")
        f.write("<div class=\"content\">\n{{ regulation_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GTEx eQTLs</button>\n")
        f.write("<div class=\"content\">\n{{ gtex_genes_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby genes</button>\n")
        f.write("<div class=\"content\">\n{{ gene_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Pubmed entries</button>\n")
        f.write("<div class=\"content\">\n{{ pubmed_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("</div>\n")

    # Variant div
    f.write("\n</div>\n\n")

    # # Genes
    f.write("<div id=\"Genes\" class=\"tabcontent2\">\n")

    # Gene buttons
    f.write("<div class=\"genetab\">\n")
    f.write("<select id=\"geneCombo\" onchange=\"showGene(this)\">\n")
    f.write("<option value=\"\">Select gene</option>\n")
    for i in range(0,len(gene_names)):
        f.write("<option value=\"%s\">%s</option>\n" %(gene_names[i],gene_names[i]))
    f.write("</select>\n")
    f.write("</div>\n\n")
    f.write("<div id=\"space\"></div>\n")

    # Gene details
    #f.write("<div class=\"genetab2\">\n")
    for i in range(0,len(gene_names)):
        f.write("<div id=\"%s\" class=\"genetabcontent\">\n" % gene_names[i])

        f.write("<button type=\"button\" class=\"collapsible\">Gene details</button>\n")
        f.write("<div class=\"content\">\n{{gene_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Gene Onthology Annotation</button>\n")
        f.write("<div class=\"content\">\n{{go_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">UniProt Annotations</button>\n")
        f.write("<div class=\"content\">\n{{uniprot_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GWAS Associations</button>\n")
        f.write("<div class=\"content\">\n{{gwas_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">eQTLs</button>\n")
        f.write("<div class=\"content\">\n{{gtexVariants_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GXA</button>\n")
        f.write("<div class=\"content\">\n{{gxa_heatmap_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Mouse Models</button>\n")
        f.write("<div class=\"content\">\n{{mouse_table_%s }}\n</div>\n" % gene_names[i])
        f.write("<div id=\"space\"></div>\n")

        f.write("</div>\n\n")

    f.write("</div>\n")

    f.write("</body>\n\n")

    f.write("<script>\n")
    f.write("function openTab2(evt, tabName) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"tabcontent2\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ntablinks = document.getElementsByClassName(\"tablinks2\");\nfor (i = 0; i < tablinks.length; i++) {\ntablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n}\ndocument.getElementById(tabName).style.display = \"block\";\nevt.currentTarget.className += \" active\";\n}\ndocument.getElementById(\"defaultTabOpen\").click();\n")

    f.write("function showGene(obj) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"genetabcontent\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ndocument.getElementById(obj.value).style.display = \"block\";\n}\n")

    f.write("function openTab(evt, tabName) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"tabcontent\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ntablinks = document.getElementsByClassName(\"tablinks\");\nfor (i = 0; i < tablinks.length; i++) {\ntablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n}\ndocument.getElementById(tabName).style.display = \"block\";\nevt.currentTarget.className += \" active\";\n}\n")

    f.write("var coll = document.getElementsByClassName(\"collapsible\");\nvar i;\n\n")
    f.write("for (i = 0; i < coll.length; i++) {\ncoll[i].addEventListener(\"click\", function() {\nthis.classList.toggle(\"active\");\nvar content = this.nextElementSibling;\nif (content.style.display === \"block\") {\ncontent.style.display = \"none\";\n} else {\ncontent.style.display = \"block\";\n}\n});\n}\n\n")

    f.write("</script>\n</html>\n")

    f.close()

# ======================================================================================================================
def generateVarTemplate_orig(mapping_names,fname):
    f = open(fname,"w")
    f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n")
    f.write("<style>\n\nth{ text-align: left; }\n")
    f.write("#common {\nborder-collapse: collapse;\nwidth: 100%;\ncolumn-width=auto;\ncolor: black;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 15px;\n}\n\n")
    f.write("div.general {\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\nmargin:  10px 10px 5px 10px\n}\n")
    f.write("h1 {\nbackground-color: #8b9bc1;\nmargin: 0;\npadding: 4px 8px 2px 24px;\n-webkit-border-radius: 8px 0 8px 0;\nline-height: 1em;\ndisplay: block;\n-webkit-margin-before: 0px;\n-webkit-margin-after: 0px;-webkit-margin-start: 0px;\n-webkit-margin-end: 0px;\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 25px;\nfont-weight: bold;\n}\n\n")
    f.write("div.big_box {\nbackground-color: white;\nwidth: 98%;\nborder: 1px solid navy;\ndisplay: block;\n\nbackground-color: #fff;\npadding: 0;\nmargin: 0;\n\n/* Margins */\nmargin: 0 auto 1em;\n\n/* Rounded edges */\nborder-bottom-right-radius: 8px;\nborder-top-left-radius: 8px;\n\n/* Shadows around boxes*/\nbox-shadow: 4px 4px 10px #BCBCCC;\n\n/* Setting fonts */\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\n}\n\n")
    f.write(".tab {\noverflow: hidden;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\n}\n")
    f.write(".tab button {\nbackground-color: inherit;\nfloat: left;\nborder: none;\noutline: none;\ncursor: pointer;\npadding: 14px 16px;\ntransition: 0.3s;\n}\n")
    f.write(".tab button:hover {\nbackground-color: #ddd;\n}\n")
    f.write(".tab button.active {\nbackground-color: #ccc;\n}\n")
    f.write(".tabcontent {\ndisplay: none;\npadding: 6px 12px;\nborder: 1px solid #ccc;\nborder-top: none;\n}\n")

    # collapsible styles
    f.write(".collapsible {\nbackground-color: #33b2ff;\ncolor: white;\ncursor: pointer;\npadding: 18px;\nwidth: 100%;\nborder: none;\ntext-align: center;\noutline: none;\nfont-size: 15px;\nborder-radius: 12px;\n}\n\n")
    f.write(".active, .collapsible:hover {\nbackground-color: #3399ff;\n}\n\n")
    f.write(".content {\npadding: 0 18px;\ndisplay: none;\noverflow: hidden;\nbackground-color: #f1f1f1;\n}\n\n")
    # ------------------

    f.write("#space\n{width: 100%;\n height: 5px;\n margin: 0px;\n padding: 0px;\n border: 0px;\nbackground: #FFFFFF;\nclear:both;\n}\n\n") 

    f.write("</style>\n\n")

    f.write("</head>\n<body bgcolor=\"#E6E6FA\" class=body>\n\n<h2>Variant mappings</h2>\n")
    f.write("<div class=\"tab\">")
    for i in range(0,len(mapping_names)):
        if i==0:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\" id=\"defaultOpen\">%s</button>" % (mapping_names[i],mapping_names[i]))
        else:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\">%s</button>" % (mapping_names[i],mapping_names[i]))
    f.write("</div>\n")
    for i in range(0,len(mapping_names)):
        f.write("<div id=\"%s\" class=\"tabcontent\">" % mapping_names[i])

        f.write("<button type=\"button\" class=\"collapsible\">Variant Details</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ variant_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby variants associated with phenotypes</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ phenotype_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby GWAS signals</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ gwas_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Predicted consequences</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ vep_table%d  }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">1kG allele frequencies</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ population_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">gnomAD allele frequencies</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ gnomad_table%d  }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">ENSEMBL Regulation</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ regulation_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GTEx eQTLs</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ gtex_genes_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Nearby genes</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ gene_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Pubmed entries</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{ pubmed_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        #f.write("<div id=\"General\" class=\"big_box\"><h1>Pubmed entries</h1>\n{{ pubmed_table%d }}\n</div>\n" %i)
        f.write("</div>\n")

    f.write("</body>\n<script>\nfunction openTab(evt, tabName) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"tabcontent\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ntablinks = document.getElementsByClassName(\"tablinks\");\nfor (i = 0; i < tablinks.length; i++) {\ntablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n}\ndocument.getElementById(tabName).style.display = \"block\";\nevt.currentTarget.className += \" active\";\n}\ndocument.getElementById(\"defaultOpen\").click();\n")

    f.write("var coll = document.getElementsByClassName(\"collapsible\");\nvar i;\n\n")
    f.write("for (i = 0; i < coll.length; i++) {\ncoll[i].addEventListener(\"click\", function() {\nthis.classList.toggle(\"active\");\nvar content = this.nextElementSibling;\nif (content.style.display === \"block\") {\ncontent.style.display = \"none\";\n} else {\ncontent.style.display = \"block\";\n}\n});\n}\n\n")

    f.write("</script>\n</html>\n")

    f.close()

# ======================================================================================================================

def generateGeneTemplate(gene_names,fname):
    f = open(fname,"w")
    f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n")
    f.write("<style>\n\nth{ text-align: left; }\n")
    f.write("#common {\nborder-collapse: collapse;\nwidth: 100%;\ncolumn-width=auto;\ncolor: black;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 15px;\n}\n\n")
    f.write("div.general {\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\nmargin:  10px 10px 5px 10px\n}\n")
    f.write("h1 {\nbackground-color: #8b9bc1;\nmargin: 0;\npadding: 4px 8px 2px 24px;\n-webkit-border-radius: 8px 0 8px 0;\nline-height: 1em;\ndisplay: block;\n-webkit-margin-before: 0px;\n-webkit-margin-after: 0px;-webkit-margin-start: 0px;\n-webkit-margin-end: 0px;\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 25px;\nfont-weight: bold;\n}\n\n")
    f.write("div.big_box {\nbackground-color: white;\nwidth: 98%;\nborder: 1px solid navy;\ndisplay: block;\n\nbackground-color: #fff;\npadding: 0;\nmargin: 0;\n\n/* Margins */\nmargin: 0 auto 1em;\n\n/* Rounded edges */\nborder-bottom-right-radius: 8px;\nborder-top-left-radius: 8px;\n\n/* Shadows around boxes*/\nbox-shadow: 4px 4px 10px #BCBCCC;\n\n/* Setting fonts */\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\n}\n\n")
    f.write("* {box-sizing: border-box}")
    f.write(".tab {\nfloat: left;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\nwidth: 10%;\nheight: 200px;\n}\n")
    f.write(".tab button {\ndisplay: block;\nbackground-color: inherit;\ncolor: black;\npadding: 22px 16px;\nwidth: 100%;\nborder: none;\noutline: none;\ntext-align: left;\ncursor: pointer;\ntransition: 0.3s;\nfont-size: 12px;\n}\n")
    f.write(".tab button:hover {\nbackground-color: #ddd;\n}\n.tab button.active {\nbackground-color: #ccc;\n}\n.tabcontent {\nfloat: left;\npadding: 0px 12px;\nborder: 1px solid #ccc;\nwidth: 85%;\nheight: 98%;\n}\n")

    # collapsible styles
    f.write(".collapsible {\nbackground-color: #33b2ff;\ncolor: white;\ncursor: pointer;\npadding: 18px;\nwidth: 100%;\nborder: none;\ntext-align: center;\noutline: none;\nfont-size: 15px;\nborder-radius: 12px;\n}\n\n")
    f.write(".active, .collapsible:hover {\nbackground-color: #3399ff;\n}\n\n")
    f.write(".content {\npadding: 0 18px;\ndisplay: none;\nalign:center\noverflow: hidden;\nbackground-color: #f1f1f1;\n}\n\n")
    # ------------------

    f.write("#space\n{width: 100%;\n height: 5px;\n margin: 0px;\n padding: 0px;\n border: 0px;\nbackground: #FFFFFF;\nclear:both;\n}\n\n")
    f.write("</style>\n")

    f.write("</head>\n<body bgcolor=\"#E6E6FA\" class=body>\n\n<h2>Genes</h2>\n")
    f.write("<div class=\"tab\">\n")
    for i in range(0,len(gene_names)):
        if i==0:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\" id=\"defaultOpen\">%s</button>\n" % (gene_names[i],gene_names[i]))
        else:
            f.write("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\">%s</button>\n" % (gene_names[i],gene_names[i]))
    f.write("</div>\n\n")
    for i in range(0,len(gene_names)):
        f.write("<div id=\"%s\" class=\"tabcontent\">\n" % gene_names[i])


        f.write("<button type=\"button\" class=\"collapsible\">Gene details</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{gene_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Gene Onthology Annotation</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{go_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">UniProt Annotations</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{uniprot_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GWAS Associations</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{gwas_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">eQTLs</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{gtexVariants_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">GXA</button>\n")
        #f.write("<div id=\"General\" class=\"content\">\n{{gxa_table%d }}\n</div>\n" %i)
        f.write("<div id=\"General\" class=\"content\">\n{{gxa_heatmap%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        f.write("<button type=\"button\" class=\"collapsible\">Mouse Models</button>\n")
        f.write("<div id=\"General\" class=\"content\">\n{{mouse_table%d }}\n</div>\n" %i)
        f.write("<div id=\"space\"></div>\n")

        #f.write("<div id=\"General\" class=\"big_box\"><h1>OMIM annotations</h1>\n<div id=\"variant\" class=\"general\" >\n{{ vep_table%d }}\n</div>\n</div>\n" %i)
        #f.write("<div id=\"General\" class=\"big_box\"><h1>Gene Expression</h1>\n{{ gnomad_table%d }}\n</div>\n" %i)

        #f.write("<div id=\"General\" class=\"big_box\"><h1>Mouse Models</h1>\n{{mouse_table%d }}\n</div>\n" %i)
        f.write("</div>\n\n")

    f.write("</body>\n<script>\nfunction openTab(evt, tabName) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"tabcontent\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ntablinks = document.getElementsByClassName(\"tablinks\");\nfor (i = 0; i < tablinks.length; i++) {\ntablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n}\ndocument.getElementById(tabName).style.display = \"block\";\nevt.currentTarget.className += \" active\";\n}\ndocument.getElementById(\"defaultOpen\").click();\n")

    f.write("var coll = document.getElementsByClassName(\"collapsible\");\nvar i;\n\n")
    f.write("for (i = 0; i < coll.length; i++) {\ncoll[i].addEventListener(\"click\", function() {\nthis.classList.toggle(\"active\");\nvar content = this.nextElementSibling;\nif (content.style.display === \"block\") {\ncontent.style.display = \"none\";\n} else {\ncontent.style.display = \"block\";\n}\n});\n}\n\n")

    f.write("</script>\n</html>\n")

    f.close()
