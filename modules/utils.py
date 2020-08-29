import logging
import re
import jinja2
import os
import sys
import pandas as pd
import tempfile as tf
import plotly.graph_objects as go

from . import config

LOGGER=logging.getLogger(__name__)
# LOGGER.setLevel(logging.DEBUG)
# ch=logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
# formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(funcName)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
# ch.setFormatter(formatter)
# LOGGER.addHandler(ch)

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
    LOGEER.debug("temp file: out.name")
    out.close()
    with open (out.name, "r") as f:
        data=f.readlines()
    f.close()
#    if os.path.isfile(out.name):
#        os.remove(out.name)

    return "".join(data)

# ======================================================================================================================

def makeLink(url,text):
    '''
    Create an HTML liink
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

# ======================================================================================================================

def chunks(L,n):
    """Yields successive n-sized chunks from L"""
    for i in range(0, len(L), n):
        yield L[i:i + n]

# ======================================================================================================================

def getLocationString(chrom,pos,ref,alt):
    '''
    Input: chrom, pos, ref, alt
    Output: VEP-style string
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
        LOGGER.error("Wrong allele encoding ref=%s, alt=%s" %(ref,alt),file=sys.stderr)
        return None

    return chrom+":"+str(start)+"-"+str(end)+":"+allele1+"/"+allele2

# ======================================================================================================================

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

# ======================================================================================================================

def generateHTML(templateFile, data):
    templateLoader = jinja2.FileSystemLoader(searchpath="/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(templateFile)
    return template.render(data)

# ======================================================================================================================

def checkFiles(files):
    '''
    Input: list of file names
    Output: True if files exist, False otherwise
    '''

    f=True
    for f in files:
        if not os.path.isfile(f):
            LOGGER.error("File %s does not exist" % f)
            f=False

    return f

# ======================================================================================================================

def createDir(path):
    if os.path.exists(path):
        #LOGGER.info("%s already exists" % path)
        return path
    try:
        os.mkdir(path)
    except OSError:
        #LOGGER.error("Creation of the directory %s failed" % path)
        return None
    else:
        return path

# ======================================================================================================================

def getGxaURL(ID):
    return config.GXA_URL_PREFIX+ID+config.GXA_URL_SUFFIX

# ======================================================================================================================

def generateTemplate(mapping_names,gene_names,fname):
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
    #f.write(".genetab {\nfloat: left;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\nwidth: 15%;\nheight: 200px;\n}\n")
#    f.write(".genetab {\noverflow: hidden;\nborder: 1px solid #ccc;\nbackground-color: #f1f1f1;\n}\n")
#    f.write(".genetab button {\nbackground-color: inherit;\nfloat: left;\nborder: none;\noutline: none;\ncursor: pointer;\npadding: 14px 16px;\ntransition: 0.3s;\n}\n")
#    f.write(".genetab button:hover {\nbackground-color: #ddd;\n}\n")
#    f.write(".genetab button.active {\nbackground-color: #ccc;\n}\n")

#    f.write(".genetab2 {\nfloat: right;\nborder: 1px solid #ccc;\nbackground-color: #0101f1;\nwidth: 80%;\n}\n")
    f.write(".tabcontent {\ndisplay: none;\npadding: 6px 12px;\nborder: 1px solid #ccc;\nborder-top: none;\n}\n")
#    f.write(".tabcontent2 {\ndisplay: none;\npadding: 6px 12px;\nborder: 1px solid #ccc;\nborder-top: none;\n}\n")
    f.write(".genetabcontent {\ndisplay: none;\nfloat: left;\npadding: 0px 12px;\nborder: 1px solid #ccc;\nwidth: 98%;\nheight: 98%;\n}\n")

    # collapsible styles
    f.write(".collapsible {\nbackground-color: #33b2ff;\ncolor: white;\ncursor: pointer;\npadding: 18px;\nwidth: 100%;\nborder: none;\ntext-align: center;\noutline: none;\nfont-size: 15px;\nborder-radius: 12px;\n}\n\n")
    f.write(".active, .collapsible:hover {\nbackground-color: #3399ff;\n}\n\n")
    f.write(".content {\nleft: 50%\npadding: 0 18px;\ndisplay: none;\ntext-align: center;\noverflow: hidden;\nbackground-color: #f1f1f1;\n}\n\n")
    # ------------------

    f.write("#space\n{width: 100%;\n height: 5px;\n margin: 0px;\n padding: 0px;\n border: 0px;\nbackground: #FFFFFF;\nclear:both;\n}\n\n") 

    f.write("</style>\n\n")

    f.write("</head>\n<body bgcolor=\"#E6E6FA\" class=body>\n\n")

    f.write("<div class=\"tab\">\n")
    f.write("<button class=\"tablinks2\" onclick=\"openTab2(event, 'Variant')\" id=\"defaultTabOpen\">Variant</button>\n")
    #f.write("<button class=\"tablinks2\" onclick=\"openTab2(event, 'Variant')\">Variant</button>\n")
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

        f.write("<button type=\"button\" class=\"collapsible\">1kG allele frequencies</button>\n")
        f.write("<div class=\"content\">\n{{ population_table%d }}\n</div>\n" %i)
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
