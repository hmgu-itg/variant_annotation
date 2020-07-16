import logging
import re
import jinja2
import os
import sys

import config

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

# ======================================================================================================================

def generateHTML(templateFile, data):
    templateLoader = jinja2.FileSystemLoader(searchpath="/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(templateFile)
    return template.render(data)

# ======================================================================================================================

def createDir(path):
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

def getGxaURL(ID):
    return config.GXA_URL_PREFIX+ID+config.GXA_URL_SUFFIX

# ======================================================================================================================

def generateVarTemplate(mapping_names,fname):
    with open(fname,"w") as f:
        sys.stdout = f

        print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n<head>\n")
        print("<style>\n\nth{ text-align: left; }\n")
        print("#common {\nborder-collapse: collapse;\nwidth: 100%;\ncolumn-width=auto;\ncolor: black;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 15px;\n}\n\n")
        print("div.general {\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\nmargin:  10px 10px 5px 10px\n}\n")
        print("h1 {\nbackground-color: #8b9bc1;\nmargin: 0;\npadding: 4px 8px 2px 24px;\n-webkit-border-radius: 8px 0 8px 0;\nline-height: 1em;\ndisplay: block;\n-webkit-margin-before: 0px;\n-webkit-margin-after: 0px;-webkit-margin-start: 0px;\n-webkit-margin-end: 0px;\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 25px;\nfont-weight: bold;\n}\n\n")
        print("div.big_box {\nbackground-color: white;\nwidth: 98%;\nborder: 1px solid navy;\ndisplay: block;\n\nbackground-color: #fff;\npadding: 0;\nmargin: 0;\n\n/* Margins */\nmargin: 0 auto 1em;\n\n/* Rounded edges */\nborder-bottom-right-radius: 8px;\nborder-top-left-radius: 8px;\n\n/* Shadows around boxes*/\nbox-shadow: 4px 4px 10px #BCBCCC;\n\n/* Setting fonts */\ncolor: navy;\nfont-family: \"Times New Roman\", Times, serif;\nfont-size: 20px;\n}\n\n")
        #print("* {box-sizing: border-box}")
        # print("/* Style the tab */\n.tab {\nfloat: left;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\nwidth: 10%;\nheight: 200px;\n}\n")
        # print("/* Style the buttons inside the tab */\n.tab button {\ndisplay: block;\nbackground-color: inherit;\ncolor: black;\npadding: 22px 16px;\nwidth: 100%;\nborder: none;\noutline: none;\ntext-align: left;\ncursor: pointer;\ntransition: 0.3s;\nfont-size: 12px;\n}\n")
        # print(".tab button:hover {\nbackground-color: #ddd;\n}\n.tab button.active {\nbackground-color: #ccc;\n}\n.tabcontent {\nfloat: left;\npadding: 0px 12px;\nborder: 1px solid #ccc;\nwidth: 85%;\nheight: 98%;\n}\n")

        print(".tab {\noverflow: hidden;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\n}")
        print(".tab button {\nbackground-color: inherit;\nfloat: left;\nborder: none;\noutline: none;\ncursor: pointer;\npadding: 14px 16px;\ntransition: 0.3s;\n}\n")
        print(".tab button:hover {\nbackground-color: #ddd;\n}\n")
        print(".tab button.active {\nbackground-color: #ccc;\n}\n")
        print(".tabcontent {\ndisplay: none;\npadding: 6px 12px;\nborder: 1px solid #ccc;\nborder-top: none;\n}\n")

        print("</style>\n")

        print("</head>\n<body bgcolor=\"#E6E6FA\" class=body>\n\n<h2>Variant mappings</h2>\n")

        print("<div class=\"tab\">")
        for i in range(0,len(mapping_names)):
            if i==0:
                print("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\" id=\"defaultOpen\">%s</button>" % (mapping_names[i],mapping_names[i]))
            else:
                print("<button class=\"tablinks\" onclick=\"openTab(event, '%s')\">%s</button>" % (mapping_names[i],mapping_names[i]))
                print("</div>")

        for i in range(0,len(mapping_names)):
            print("<div id=\"%s\" class=\"tabcontent\">" % mapping_names[i])
            print("<div id=\"General\" class=\"big_box\"><h1>Variant details</h1>\n<div id=\"variant\" class=\"general\" >General information:<br>\n{{ variant_table%d }}\n</div>\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Nearby variants associated with phenotypes</h1>\n{{ phenotype_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Nearby GWAS signals</h1>\n{{ gwas_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Predicted consequences</h1>\n<div id=\"variant\" class=\"general\" >\n{{ vep_table%d }}\n</div>\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>1kG allele frequencies</h1>\n{{ population_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>ExAC allele frequencies</h1>\n{{ exac_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Regulation</h1>\n{{ regulation_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>GTEx eQTLs</h1>\n{{ gtex_gene_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Nearby genes</h1>\n{{ gene_table%d }}\n</div>\n" %i)
            print("<div id=\"General\" class=\"big_box\"><h1>Pubmed entries</h1>\n{{ pubmed_table%d }}\n</div>\n" %i)
            print("</div>\n")
        
            print("</body>\n<script>\nfunction openTab(evt, tabName) {\nvar i, tabcontent, tablinks;\ntabcontent = document.getElementsByClassName(\"tabcontent\");\nfor (i = 0; i < tabcontent.length; i++) {\ntabcontent[i].style.display = \"none\";\n}\ntablinks = document.getElementsByClassName(\"tablinks\");\nfor (i = 0; i < tablinks.length; i++) {\ntablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n}\ndocument.getElementById(tabName).style.display = \"block\";\nevt.currentTarget.className += \" active\";\n}\ndocument.getElementById(\"defaultOpen\").click();\n</script>\n</html>\n")

