'''
This collection of function was designed to create html formatted tables and
paragraphs for the variant annotation tool I've been developing.

Package contents:
* draw_gene_info


Last modified: 2015.10.01

'''

# Importing standard libraries:
import sys, os, re
import matplotlib.pyplot as plt
import pandas as pd

# Importing library for creating word cloud:
sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Django-1.8.6-py2.7.egg')
sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Jinja-1.2-py2.7-linux-x86_64.egg')
sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/wordcloud-1.2-py2.7-linux-x86_64.egg')
sys.path.append('/nfs/team144/ds26/anaconda2/lib/python2.7/site-packages')


# These functions are important for creating the gene onthology word cloud:
from wordcloud import WordCloud

# Generating a html formatted table of the gene annotation.
def draw_gene_info(gene_data):
    '''
    This function creates a html formatted table with
    general information of the given gene.

    Input: Dictionary downloaded from Ensembl.
    output: html formatted table in a string.
    '''
    # Initializing the table:
    table = "<table class=\"gene_general\" width=600px >\n"

    # This small drawing function adds each rows to the table:
    def add_row(col1, col2):
        return "\t<tr><td class=\"gene_cat\">%s:</td><td>%s</td></tr>\n" %(col1, col2)

    coordinates = "chr%s:%s-%s" %(gene_data["chromosome"], gene_data["start"], gene_data["end"])
    linked_ID = "<a href=\"http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s\">%s</a>" %(gene_data["id"], gene_data["id"])

    # Adding relevant rows:
    table += add_row("Gene name", gene_data["name"])
    table += add_row("Description", gene_data["description"])
    table += add_row("Ensembl ID", linked_ID)
    table += add_row("Genomic coordinates", coordinates)
    table += add_row("Strand", gene_data["strand"])
    table += add_row("Gene type", gene_data["type"])

    # Closing table:
    table += "</table>"

    return table

# Generating word cloud of the GO terms. + html formatted list of go terms.
def draw_GO_terms(gene_xrefs, Gene_ID):
    '''
    This function takes the gene onthology cross-references downloaded from Ensembl
    and generates a formatted list with links pointing to the GO website.

    This function also generates a wordcloud image of the used GO terms. And the
    saved image is inserted into the html
    '''

    # File name of the resulting png file. It should be fixed soon:
    filename = "%s_GO_wordcloud.png" % Gene_ID

    # Generating a list of unique GO terms (go annotation is provided for all transcripts)
    unique_GO_terms = {}
    for counter, go_term in enumerate(gene_xrefs["GO"]):
            unique_GO_terms[go_term[0]] = go_term[1]

    # Extracting individual words from all unique go terms:
    words = []
    for counter, go_term in enumerate(unique_GO_terms.keys()):
        words.extend(go_term.replace("-"," ").split(" "))

    # There are some very frequently used words that we remove:
    removed_words = ["to", "response", "binding", "protein", "of"]
    words = [i for i in words if i not in removed_words]

    # Let's generate frequencies:
    freq = {}
    for word in words:
        try:
            freq[word] += 1
        except:
            freq[word] = 1

    # number of different counts for gradient:
    uniq_count = len(list(set(freq.values())))

    # Return an rgb code for each counts keeping in mind that even the
    # rarest group should be readable.
    colors = dict(zip(list(set(freq.values())), linear_gradient(start_hex = '#8b9bc1',
                    finish_hex = "#a9b5d1",
                    n = uniq_count)))

    # Now we have to update out dictionary with the colors:
    words_with_colors = {}

    # Looping through all words and assigning a color based on their
    # frequencies:
    for word in words:
        words_with_colors[word] = colors[freq[word]]

    # Calling globalizer function:
    gobal_dict(words_with_colors)

    # Initializing wordcloud object:
    wc = WordCloud(background_color="white",
                   max_words=2000,
                   color_func=get_color,
                   prefer_horizontal = 0.5)

    html = ""
    if len(words) > 0:
        # generate word cloud:
        wc.generate(" ".join(words))

        # Saving figure:
        wc.to_file("./genes/"+filename)

        # Second let's creating a list of the Gene onthology terms with links:
        html += "<div id=\"GOfig\" style=\"text-align: center\" ><img src=\"%s\" alt=\"Go-term word cloud\" style=\"width:304px;height:228px\"></div></br>\n" % filename

    html += "<div id=\"full_GO\" class=\"full_GO\">\n"
    for counter, go_term in enumerate(unique_GO_terms.keys()):
        html += "\t<a href=\"http://amigo.geneontology.org/amigo/term/%s\">%s</a>,\n" %(unique_GO_terms[go_term], go_term)

    html += "</div>"
    return html

# Formatting OMIM annotations:
def draw_OMIM_annotation(OMIM_annotation):
    '''
    This function creates a html formatted text from the OMIM annotation.

    Input: returned OMIM annotation

    return: properly formatted html string.
    '''
    html = ""

    # If the gene did not have any MIM entry, we do accordingly:
    ## Check if the returned value is a dictionary or a string:
    if type(OMIM_annotation) == str or type(OMIM_annotation) == unicode:
        html = '<div class="missing_data">%s</div>' % OMIM_annotation
        return html

    for ID in OMIM_annotation.keys():
        # Generating header:
        header = "%s <a href=\"http://www.omim.org/entry/%s\">Link</a>" % (OMIM_annotation[ID]["title"], ID)

        # Generating the body of the OMIM entry:
        #print OMIM_annotation[ID]

        # Processing text:
        text = ""
        for field in OMIM_annotation[ID]["text"].keys():
            text += "\t\t<div class=\"OMIM_cat\">%s:</div><div class=\"OMIM_desc\">\n" % field
            text += "\t\t\t<span>%s<br></span>\n" %  OMIM_annotation[ID]["text"][field].replace("\n", "<br>")
            text += "\t\t</div>\n"

        # Processing variations:
        if len(OMIM_annotation[ID]["variations"]) > 1:
            variation_html = "\t\t<div class=\"OMIM_cat\">Variations:</div><div class=\"OMIM_desc\">\n"
            for var in OMIM_annotation[ID]["variations"].keys():
                variation_html += "\t\t\t<span><b>%s: </b></span><span>%s<br></span>\n" %(var, OMIM_annotation[ID]["variations"][var])
            variation_html += "\t\t</div>\n"

            text += variation_html

        # Combining header and body together:
        html += "<div class=\"OMIM\" onclick=\"toggle()\">%s\n\t<div style=\"display: none;\">\n%s\n\t</div>\n</div>" %(header, text)

    return html

# helper function to generate color gradient:
def hex_to_RGB(hex):
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i+2], 16) for i in range(1,6,2)]

# helper function to generate color gradient:
def linear_gradient(start_hex, finish_hex="#FFFFFF", n=307):
    ''' returns a gradient list of (n) colors between
        two hex colors. start_hex and finish_hex
        should be the full six-digit color string,
        inlcuding the number sign ("#FFFFFF") '''

    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)

    # Initilize a list of the output colors with the starting color
    RGB_list = ["rgb(%s, %s, %s)" % (s[0], s[1], s[2])]

    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append("rgb(%s, %s, %s)" % (curr_vector[0], curr_vector[1], curr_vector[2]))

    RGB_list.reverse()
    return RGB_list

# A function to generate global variable:
def gobal_dict(words):
    global word_dict
    word_dict = words

# This is a helper function called by the wordcloud script.
# Uses the global variable wods_with_colos dictionary.
def get_color(word=None, font_size=None, position=None,
                      orientation=None, font_path=None, random_state=None):
    try:
        rgb_string = word_dict[word]
    except:
        rgb_string = "rgb(255, 0, 0)"

    return rgb_string

# Adding Pubmed references into uniprot text
def add_pubmed_ref(match):
    pubmed_IDs = re.findall("PubMed:(\d+)", match.group(0), flags=0)
    return_value = ""
    for r in pubmed_IDs:
        return_value += "<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/%s\">Ref</a> " % r
    return return_value

# Adding OMIM references into uniprot text
def add_OMIM_ref(match):
    pubmed_IDs = re.findall("MIM:(\d+)", match.group(0), flags=0)
    return_value = ""
    for r in pubmed_IDs:
        return_value += "<a href=\"http://www.omim.org/entry/%s\">Ref</a> " % r
    return return_value

# A function to prepare a nice html formatted text from the Uniprot annotation
def draw_Uniprot_table(Uniprot_data):
    '''
    This function was written to format data retrieved from Uniprot.
    Each entry is linked to uniprot in the header. Non empty fields
    are displayed.

    To do: format references... OMIM/Pubmed. We have to use regular
    expressoins...

    '''
    # Generating html formatted document based on the downloaded uniprot data:
    table = ""

    # Looping through all uniprot entries associated with the given gene:
    for entry in Uniprot_data.keys():
        # Generate div header:
        URL = "%s (<a href=\"http://www.uniprot.org/uniprot/%s\">%s</a>)" %(entry,
                               Uniprot_data[entry]["Entry"], Uniprot_data[entry]["Name"])
        table += "<div class=\"Uniprot_entry\" onclick=\"toggle()\" id=\"%s\"> %s</br>\n" % (Uniprot_data[entry]["Name"], URL)
        table += "<div class=\"Uniprot_hidden\" id=\"%s\" style=\"display: none;\">\n" % Uniprot_data[entry]["Name"]
        # formatting text:
        for field in Uniprot_data[entry].keys():
            # We are excluding the basic fields and empty categories:
            if field == "Name" or field == "Entry" or Uniprot_data[entry][field] == "":
                continue

            # Removing the first part:
            description = " ".join(Uniprot_data[entry][field].split(": ")[1:])

            # Replacing references with links to OMIM or pubmed:
            description = re.sub("\{.+?\}", add_pubmed_ref, description)
            description = re.sub("\[MIM:\d+\]", add_OMIM_ref, description)

            # Creating html formatted text:
            table += "\t<div class=\"UnipCategory\">%s:</div><div class=\"UnipDescription\">%s</div>\n" %(field, description)


        table += "</div></div>\n"
    table += "<br>"

    return table

# Draw gwas hits associated with the gene:
def draw_gwas_gene_table(GWAS_data):
    '''
    Creates a html formatted table of the gwas hits associated with the gene.

    Input: array with gwas hits, generated by the dedicated function.
    Output: html formatted string
    '''

    ## Check if the returned value is a dictionary or a string:
    if type(GWAS_data) == str or type(GWAS_data) == unicode:
        html = '<div class="missing_data">%s</div>' % GWAS_data
        return html

    # Initializing GWAS table:
    html = '<table id="gene_GWAS_table" class="gene_GWAS_table">\n'
    html += '\t<tr class="table_header"><td>rsID</td><td>Consequence</td><td>Reported gene</td><td>Trait</td><td>Publication</td></tr>\n'

    # Looping through all associated signals:
    for signal in GWAS_data:
        html += "\t<tr>"

        # Adding rsID:
        try:
            html += '<td><a href="https://www.ebi.ac.uk/gwas/search?query=%s">%s</a></td>' %(signal["rsID"], signal["rsID"])
        except:
            html += '<td>%s</td>' % signal["rsID"]

        # Adding consequence:
        try:
            html += '<td>%s</td>' % signal["Consequence"]
        except:
            html += '<td>NA</td>'

        # Adding reported gene:
        try:
            html += '<td>%s</td>' % signal["Reported genes"]
        except:
            html += '<td>NA</td>'

        # Adding trait:
        try:
            html += '<td>%s</td>' % signal["trait"]
        except:
            html += '<td>NA</td>'

        # Adding publication:
        try:
            html += '<td><a href="%s">%s</a></td>' % (signal["URL"], signal["Author"])
        except:
            html += '<td>NA</td>'

        html += "</tr>\n"
    html += "</table>\n"
    return html

# Creating a html formatted table of the gene expression atlas data:
def draw_GXA_table(GXA_levels, general):
    '''
    This function creates a html formatted table of the GXA data.
    If there was no data, it will print that out.

    It also creates a heatmap of the expression levels pooled togeteher
    for tissue types.
    '''

    if type(GXA_levels) == str:
        return "No gene expression informatin at Gene Expression Atlas<br>\n"

    # Initializing the top row of the field:
    table = 'Top 10 highest expressing tissues from the most important sources in <a href="http://www.ebi.ac.uk/gxa/genes/%s">Gene Expression Atlas (EBI)</a>:<br><br>\n' % general["id"]

    # Looping through all keys of the table:
    for index, GXA_list in enumerate(GXA_levels.values()):
        table += '<div class="GXA_levels" onclick="toggle()" style="margin-left:10px; margin-bottom:2px">%s:\n' % GXA_levels.keys()[index]
        table += '<div class="GXA_hidden" style="display: none;">\n'
        table += '<table class="gene_GXA_table" style="width: 400px">\n'
        table += '\t<tr class="table_header"><td>Tissue</td><td>Expression level</td></tr>\n'

        for pair in GXA_list:
            #print "%s %s" %(GXA_dict[value], value)
            table += '\t<tr><td>%s</td><td>%s</td></tr>\n' % (pair[0], pair[1])

        table += '</table>\n</div></div>\n'

    return table

def draw_GXA_heatmap(GXA_data):
    '''
    Based on the returned Gene Expression Atlas data, this function
    creates a heatmap (svg format) of the expression values for each
    tissue type.

    The returned string is a html formatted string with the heatmap
    '''

    # If there
    if type(GXA_data) == str:
        return ""

    # The following dictionaries will be used to assign tissue category and colors for each datapoint:
    tissueMap = {
        'EBV-transformed lymphocyte' : "Hematopoietic system",
        'adipose' : "Connective tissue",
        'adipose tissue' : "Connective tissue",
        'adrenal' : "Gland",
        'adrenal gland' : "Gland",
        'amygdala' : "Nervous system",
        'animal ovary' : "Gland",
        'anterior cingulate cortex (BA24) of brain' : "Nervous system",
        'aorta' : "Cardiovascular system",
        'appendix' : "Digestive system",
        'arm muscle' : "Musculoskeletal system",
        'artery' : "Cardiovascular system",
        'atrial appendage of heart' : "Cardiovascular system",
        'bladder' : "Urogenital system",
        'bone marrow' : "Hematopoietic system",
        'brain' : "Nervous system",
        'breast' : "Gland",
        'breast (mammary tissue)' : "Gland",
        'caudate (basal ganglia) of brain' : "Nervous system",
        'caudate nucleus' : "Nervous system",
        'cerebellar hemisphere of brain' : "Nervous system",
        'cerebellum' : "Nervous system",
        'cerebral cortex' : "Nervous system",
        'cerebral meninges' : "Nervous system",
        'cervix' : "Urogenital system",
        'colon' : "Digestive system",
        'coronary artery' : "Cardiovascular system",
        'cortex of kidney' : "Gland",
        'diaphragm' : "Musculoskeletal system",
        'diencephalon' : "Nervous system",
        'duodenum' : "Digestive system",
        'dura mater' : "Nervous system",
        'ectocervix' : "Urogenital system",
        'endometrium' : "Urogenital system",
        'epididymis' : "Urogenital system",
        'esophagus' : "Digestive system",
        'esophagus muscularis mucosa' : "Digestive system",
        'eye' : "Sense organ",
        'fallopian tube' : "Urogenital system",
        'frontal cortex (BA9)' : "Nervous system",
        'frontal lobe' : "Nervous system",
        'gall bladder' : "Urogenital system",
        'gastroesophageal junction' : "Digestive system",
        'globus pallidus' : "Nervous system",
        'heart' : "Cardiovascular system",
        'heart left ventricle' : "Cardiovascular system",
        'hippocampus' : "Nervous system",
        'hypothalamus' : "Nervous system",
        'kidney' : "Gland",
        'large intestine' : "Digestive system",
        'left atrium' : "Cardiovascular system",
        'left kidney' : "Gland",
        'left renal cortex' : "Gland",
        'left renal pelvis' : "Gland",
        'left ventricle' : "Cardiovascular system",
        'leg muscle' : "Musculoskeletal system",
        'leukemia cell line' : "Hematopoietic system",
        'leukocyte' : "Hematopoietic system",
        'liver' : "Gland",
        'locus coeruleus' : "Nervous system",
        'lung' : "Respiratory system",
        'lymph node' : "Gland",
        'medulla oblongata' : "Nervous system",
        'middle frontal gyrus' : "Nervous system",
        'middle temporal gyrus' : "Nervous system",
        'minor salivary gland' : "Gland",
        'mitral valve' : "Cardiovascular system",
        'mucosa of esophagus' : "Digestive system",
        'nucleus accumbens (basal ganglia)' : "Nervous system",
        'occipital cortex' : "Nervous system",
        'occipital lobe' : "Nervous system",
        'olfactory apparatus' : "Sense organ",
        'ovary' : "Urogenital system",
        'pancreas' : "Digestive system",
        'parietal lobe' : "Nervous system",
        'parotid gland' : "Gland",
        'penis' : "Urogenital system",
        'pineal gland' : "Nervous system",
        'pituitary gland' : "Gland",
        'placenta' : "Cardiovascular system",
        'prefrontal cortex' : "Nervous system",
        'prostate' : "Gland",
        'prostate gland' : "Gland",
        'pulmonary valve' : "Cardiovascular system",
        'putamen' : "Nervous system",
        'putamen (basal ganglia)' : "Nervous system",
        'rectum' : "Digestive system",
        'renal cortex' : "Gland",
        'renal pelvis' : "Gland",
        'right renal cortex' : "Gland",
        'right renal pelvis' : "Gland",
        'salivary gland' : "Gland",
        'seminal vesicle' : "Urogenital system",
        'sigmoid colon' : "Digestive system",
        'skeletal muscle' : "Musculoskeletal system",
        'skin' : "Skin",
        'skin of lower leg' : "Skin",
        'skin of suprapubic region' : "Skin",
        'small intestine' : "Digestive system",
        'smooth muscle' : "Musculoskeletal system",
        'spinal cord' : "Nervous system",
        'spinal cord (cervical c-1)' : "Nervous system",
        'spleen' : "Hematopoietic system",
        'stomach' : "Digestive system",
        'subcutaneous adipose tissue' : "Connective tissue",
        'submandibular gland' : "Gland",
        'substantia nigra' : "Nervous system",
        'temporal lobe' : "Nervous system",
        'terminal ileum of small intestine' : "Digestive system",
        'testis' : "Urogenital system",
        'thalamus' : "Nervous system",
        'throat' : "Respiratory system",
        'thymus' : "Gland",
        'thyroid' : "Gland",
        'tibial artery' : "Cardiovascular system",
        'tibial nerve' : "Nervous system",
        'tongue' : "Digestive system",
        'tonsil' : "Hematopoietic system",
        'trachea' : "Respiratory system",
        'transformed fibroblast' : "Connective tissue",
        'transverse colon' : "Digestive system",
        'triscuspid valve' : "Cardiovascular system",
        'trunk muscle' : "Musculoskeletal system",
        'umbilical cord' : "Cardiovascular system",
        'urinary bladder' : "Urogenital system",
        'uterus' : "Urogenital system",
        'vagina' : "Urogenital system",
        'vas deferens' : "Urogenital system",
        'visceral (omentum) adipose tissue' : "Connective tissue",
        'whole blood' : "Hematopoietic system"
    }

    color_gradient = {
        1 : "#8B9BC1",
        0.8 : "#97A6C7",
        0.7 : "#A4B1CE",
        0.6 : "#B1BCD5",
        0.5 : "#BEC7DC",
        0.4 : "#CBD2E3",
        0.3 : "#D8DDEA",
        0.2 : "#E5E8F1",
        0.1 : "#F2F3F8",
        0.0 : "#FFFFFF",
        "NaN" : "#E0E0E0"}

    # Combining gxa data together based on tissue category:
    tissues = []
    for tgroup in set(tissueMap.values()):
        tissues = [tissue for tissue, group in tissueMap.items() if group == tgroup and tissue in GXA_data.columns]
        GXA_data[tgroup] = GXA_data[tissues].max(axis=1, skipna=True, level=None, numeric_only=True,)

    # Combining data together:
    experiments = GXA_data["Experiment"]
    toHeatMap = GXA_data[list(set(tissueMap.values()))]
    toHeatMap.index = GXA_data["Experiment"].tolist()

    # Creating the svg formatted heatmap:
    x = 200 # x coordinate of the first column
    y = 120 # y coordinate of the first row
    size = 30 # Size of the block
    rotation = -45 # Degree of rotation of the tissue type

    # Initializing the svg block:
    heatmap = '<div style="margin-top: 25px">Heatmap summarizing relative expression levels in tissue categories:<br><svg width="650" height="%s">\n' %(130 + len(toHeatMap)*30)

    # Printing tissue categories:
    for tissue in set(tissueMap.values()):
        x += 30
        tissue_group = ",&#013;".join([t for t, group in tissueMap.items() if group == tissue])
        heatmap += '\t<text x="%s" y="%s" fill="black" style="font-size: 15px;" transform="rotate(%s %s, %s)"><title>%s</title>%s</text>\n' % (x+30, y, rotation, x+30 ,y, tissue_group, tissue)

    # Reset coordinates:
    x = 215
    y = 130
    for experiment, row in toHeatMap.iterrows():
        # Adding row label:
        heatmap += '\t<text x="%s" y="%s" fill="black" style="font-size: 20px;">%s</text>\n' % (15, y+20, experiment)

        rowmax = toHeatMap.ix[experiment].max()
        for index, cell in enumerate(row):
            x += size
            tissue = row.index[index]
            try:
                color = color_gradient[round(cell / rowmax, 1)]
            except:
                color = color_gradient["NaN"]

            heatmap += '\t<rect x="%s" y="%s" width="%s" height="%s" style="fill:%s;stroke-width:0;fill-opacity:1" />\n' % (x, y, size, size, color)

        y += size
        x = 215

    heatmap += '</svg>\n<div style="font-size:16px;margin-left:20px;margin-right:20px"><br>'
    heatmap += '(Heatmap indicates relative expression levels for each source, where each tissue type is represented by the highest expression level.'
    heatmap += 'Move cursor over tissue category names to see which tissues are pooled in that category. '
    heatmap += 'Grey squares indicate missing values.)</div></div>\n'
    return heatmap

def draw_GTEx_variations(GTEx_data):
    '''
    This function generates a html formatted report based on the returned GTEx data.
    Input: dictionary generateg by get_GTEx_variations
    Output: html formatted string.
    '''

    ## Check if the returned value is a dictionary or a string:
    if type(GTEx_data) == str or type(GTEx_data) == unicode:
        html = '<div class="missing_data">%s</div>' % GTEx_data
        return html

    # template strings for links:
    variation_link = '<a href="http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=%s;vdb=variation">%s</a>'
    gene_link = '<a href="http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s">%s</a>'
    GTEX_link = '(<a href="http://www.gtexportal.org/home/gene/%s">GTEx site</a>)'

    # If a dictionary, initialize table:
    table = '<table class="GTEX_eQTL">\n'
    table += '\t<tr class="Header"><td>rsID</td><td>Tissue</td><td>p-value</td><td>Distance from TSS</td></tr>\n'

    flag = 0 # Flag used for coloring the background of the table
    rsID_list = {} # Keeping track how many rsIDs effect the expression of the given gene
    tissue_list = {} # Keeping track how many tissues are the expression effected by the given rsID
    gene_name = ""
    gene_ID = ""
    for rsID in GTEx_data:

        # extracting gene name and ID:
        if not gene_name:
            gene_name = GTEx_data[rsID][0]["gene_name"]
            gene_ID = GTEx_data[rsID][0]["gene_ID"]

        # in each rsID we change the background color of the table:
        color = "#E5E5F1"
        if flag == 0:
            color = "#FFFFFF"
            flag = 1
        else:
            flag = 0

        rsID_link = variation_link %(rsID, rsID)
        table += '\t<tr style="background-color: %s"><td valign="top" rowspan="%s">%s</td></tr>\n' % (color, len (GTEx_data[rsID]) + 1, rsID_link)

        rsID_list[rsID] = 1 # updating hash

        for assoc in GTEx_data[rsID]:
            # Formatting distance:
            dist = int(assoc['distance'])
            if abs(dist) > 1000:
                dist = str(dist/1000)+"kbp"
            else:
                dist = str(dist)+"bp"

            # Formatting p-value:
            (num, exp) = assoc['pvalue'].split("e")
            num = '%.3f' % round(float(num), 3) +"e"+exp

            table += '\t\t<tr style="background-color: %s"><td>%s</td><td>%s</td><td align="center">%s</td></tr>\n' % (color, assoc['tissue'], num, dist)
            tissue_list[assoc['tissue']] = 1 # Updating tissue counter

    table += '</table>\n'

    # Generate a report about the finding:
    gene_link = gene_link % (gene_ID, gene_name)
    Report = "The expression of %s was affected by %s variations in %s tissues." %(gene_link, len(rsID_list), len(tissue_list))
    Report += GTEX_link % (gene_ID)

    # Adding report and table together in a html formatted string:
    html = '<div class="GTEx variation" onclick="toggle()" style="margin-left:10px; margin-bottom:2px">%s\n' % Report
    html += '<div class="GXA_hidden" style="display: none;">\n%s\n</div></div>\n' %table

    # Return:
    return html

def draw_mouse_phenotypes(mouse_df):
    '''
    Input: dataframe if mouse phenotypes are known.
        string if returning mouse phenotype failed somewhere
    Output: html formatted string.
    '''

    ## Check if the returned value is a string or a dataframe:
    if type(mouse_df) == str or type(mouse_df) == unicode:
        html = '<div class="missing_data">%s</div>' % mouse_df
        return html

    # helper function that will be used later:
    def _get_table_line (row):
        # Processing lines:
        URL_allele = "http://www.informatics.jax.org/allele/%s" % row.Allele_ID
        try:
            phenotypes = row.Phenotypes.replace(" | ", "<br>")
        except:
            phenotypes = row.Phenotypes

        # generate line
        return  '\t<tr><td><a href="%s">%s</a></td><td>%s</td><td>%s</td><td>%s</td></tr>' % (
            URL_allele, row.Allele_name, row.Allele_type, phenotypes, row.Human_disease)

    # Generating a smaller table witht the mouse gene name, description, ID and MGI ID.
    mouse_df = mouse_df.set_index("mouse_gene_ID")
    mouse_gene_data = pd.DataFrame({"name": mouse_df.mouse_gene_name.unique(),
                                    "desc" : mouse_df.mouse_gene_description.unique(),
                                    "MGI" : mouse_df.MGI_ID.unique()},
                                  index=mouse_df.index.unique(),)

    # The table header is constant:
    table_header = '\n<table class="mouse_phenotype">\n\t<tr class="mouse_pheno_header">'\
                '<td style="width:50%">Allele name</td>'\
                '<td style="width:10%">Allele type</td>'\
                '<td style="width:15%">Phenotype</td>'\
                '<td style="width:15%">Human disease</td></tr>\n'
    html = "" # Initializing html table.

    # Now looping through all the different cases and generage tables:
    for ID in mouse_gene_data.index:
        # Printing section header:
        html += '<div id="mouse_genes" class="mouse_gene" style="font-size:17px">'\
                '<a href="http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s">%s</a>'\
                '(<a href="http://www.informatics.jax.org/marker/%s">MGI link</a>): %s<br></div>\n' % (
            ID, mouse_gene_data.loc[ID]["name"], mouse_gene_data.loc[ID]["MGI"], mouse_gene_data.loc[ID]["desc"])

        # Printing table:
        html += table_header # header
        html += "\n".join(mouse_df[mouse_df.index == ID].apply(_get_table_line, axis=1)) # rows
        html += '\n</tbody></table>\n\n' # Closing table

    return html