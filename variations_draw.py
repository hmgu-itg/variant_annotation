'''
These functions were collected to draw html formatted tables based on variation data

version: v.2.1.0 Last modified: 2015.11.23

Package includes:

'''
# Libraries for configuration:
import config

# Shared functions:
from shared import *

## Loading sys is required to load the jinja lib:
#import sys
#
## Loading this module for debug purposes:
#from pprint import pprint
#
## Importing the templating library that creates html files:
#sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Django-1.8.6-py2.7.egg')
#sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/Jinja-1.2-py2.7-linux-x86_64.egg')
#sys.path.append('/nfs/users/nfs_d/ds26/local/lib/python2.7/site-packages/wordcloud-1.2-py2.7-linux-x86_64.egg')
#sys.path.append('/nfs/team144/ds26/anaconda2/lib/python2.7/site-packages')
#
#
## Load the jinja library's namespace into the current module.
#import jinja2
#
#
##
#def draw_html(templateFile, data):
#    '''
#    This is the html templating script. It needs the template html file, and the printed data set (in dictionary format).
#
#    Output is the html string
#    '''
#
#    # search path for template files:
#    templateLoader = jinja2.FileSystemLoader( searchpath="/" )
#
#    # An environment provides the data necessary to read and
#    templateEnv = jinja2.Environment( loader=templateLoader )
#
#    # Read the template file using the environment object.
#    # This also constructs our Template object.
#    template = templateEnv.get_template( templateFile )
#
#    # Specify any input variables to the template as a dictionary.
#    templateVars = data
#
#    # Finally, process the template to produce our final text.
#    outputText = template.render( templateVars )
#
#    return outputText

#
def draw_variation_table(variation_detail):
    '''
    This function prints out a html formatted table of the allele data. It could contain various details dependeing
    on the available information.

    Anchors for Ensembl and other sources are read from config.py
    '''

    # Initializing the table:
    table = '<table class="variation" width=600px >\n'

    # This small drawing function adds each rows to the table:
    def add_row(col1, col2):
        return "\t<tr><td class=\"var_cat\">%s:</td><td>%s</td></tr>\n" %(col1, col2)

    # Adding rsID and position rows. They will contain links to Ensembl:
    if "rs" in variation_detail["rsID"]:
        rsID_row = config.ENSEMBL_VAR % (variation_detail["rsID"], variation_detail["rsID"])
    else:
        rsID_row = variation_detail["rsID"]

    position_row = config.ENSEMBL_REGION % (variation_detail["location"], variation_detail["location"])

    # Adding values to table:
    table += add_row("rsID", rsID_row)
    table += add_row("Chromosome:start-end", position_row)
    table += add_row("Allele string", variation_detail["allele_string"])
    table += add_row("MAF", variation_detail["MAF"])
    table += add_row("Consequence", variation_detail["consequence"])
    table += add_row("Variation type", variation_detail["varClass"])

    # The following fields are optinal, might not be filled:
    if "gerp" in variation_detail:
        table += add_row("gerp (average gerp)", "%s (%s)" %(variation_detail["gerp"], variation_detail["avg_gerp"]))
    if "gwava" in variation_detail:
        table += add_row("GWAVA score", variation_detail["gwava"])
    if "ancestralAllele" in variation_detail:
        table += add_row("Ancestral allele", variation_detail["ancestralAllele"])
    if "minorAllele" in variation_detail:
        table += add_row("Minor allele", variation_detail["minorAllele"])
    if "synonyms" in variation_detail and len(variation_detail["synonyms"]) > 0:
        table += add_row("Synoyms", ", ".join(variation_detail["synonyms"]))

    # Some cases, if the variation overlaps with protein coding regions, sift and polyphen scores are also available:
    if len(variation_detail["sift"]) > 0:
        try:
            table += add_row( "SIFT (score)", "%s (%s)" % (variation_detail["sift"][1],variation_detail["sift"][0] ))
        except:
            table += add_row( "SIFT (score)", "%s" % (variation_detail["sift"][0] ))
    if len (variation_detail["polyphen"]) > 0:
        try:
            table += add_row( "Polyphen (score)", "%s (%s)" % (variation_detail["polyphen"][1],variation_detail["polyphen"][0] ))
        except:
            table += add_row( "Polyphen (score)", "%s" % (variation_detail["polyphen"][0] ))

    # Phenotypes and clinical data is raraly available:
    if "phenotypes" in variation_detail.keys() and len(variation_detail["phenotypes"]) > 0:
        # For generating links for phenotype entries:
        line = []
        for phenotypes in variation_detail["phenotypes"]:
            if phenotypes["source"]  == "ClinVar":
                line.append(config.CLINVAR_LINK % (variation_detail["rsID"], phenotypes["trait"]))
            elif phenotypes["source"]  == "OMIM":
                line.append(config.OMIM_LINK % (variation_detail["rsID"], phenotypes["trait"]))
            elif phenotypes["source"]  == "NHGRI-EBI GWAS catalog":
                line.append(config.GWAS_CAT_LINK % (variation_detail["rsID"], phenotypes["trait"]))
            else:
                line.append("%s" % (phenotypes["trait"]))

        table += add_row("Phenotypes<br>(risk allele)", ",<br>".join(line))

    # If clinical significance is given for the variation we are adding extra row in the table:
    if "clinical_significance" in variation_detail.keys() and len(variation_detail["clinical_significance"]) > 0:
        table += add_row("Clinical significance ", ", ".join(variation_detail["clinical_significance"]))

    # If the variation is missense, and the codon change is given, we are adding extra row in the table:
    if "codon" in variation_detail.keys() and len(variation_detail["codon"]) > 0:
        table += add_row("Codon change", variation_detail["codon"])
    if "aminoacid" in variation_detail.keys() and len(variation_detail["aminoacid"]) > 0:
        table += add_row("Amino acid change", variation_detail["aminoacid"])


    table += "</table>"
    return table

#
def draw_gwas_position_table(gwas_table):
    '''
    This function creates a html formatted string of all gwas hits generated by get_gwas_hits_position() function.
    Links to the GWAS catalog and to the Pubmed entries are stored in the config.py file.
    '''

    # Check if the returned value is a dictionary or a string:
    if CheckReturnedValue(gwas_table): return CheckReturnedValue(gwas_table)

    # Drawing header:
    table = "<table class=\"gwas_table\">\n\t<tr><td class=\"table_header\">rsID</td><td class=\"table_header\">SNPID</td><td class=\"table_header\">Distance</td><td class=\"table_header\" width=\"35%\">Trait</td><td class=\"table_header\">p-value</td><td class=\"table_header\">PMID</td></tr>\n"

    # Blank row to be filled up with data:
    table_row = "\t<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n"

    for hit in gwas_table:
        # Filling up fields:
        if hit['PMID'] == "NA":
            table += table_row % (hit['rsID'],hit['SNPID'], hit['distance'],hit['trait'],hit['p-value'],hit['PMID'])
        else:
            table += table_row % (config.GWAS_CAT_LINK %
                                    (hit['rsID'], hit['rsID']),
                                    hit['SNPID'],
                                    hit['distance'],
                                    hit['trait'],
                                    hit['p-value'],
                                    config.PUBMED_LINK %(hit['PMID'],hit['PMID'])
                                )

    table += "</table>"

    return table

# Drawing function to generate a table with all 1000Genome frequencies:
def draw_freq_table(pops, allele_string):

    table ='<br>\n<div id="1000_genomes" class="general" > Phase 3 frequencies from the 1000 Genomes project:<br></div>\n'

    # If the returned value is a string we report the returned value:
    if CheckReturnedValue(pops): return CheckReturnedValue(pops)


    # Find alleles at first:
    alleles = allele_string.split("/")

    # Populations of the 1000 Genomes, grouped by continents:
    population_groups = {
            "AFR" : ["ACB", "ASW", "ESN", "LWK", "MAG", "MSL", "YRI"],
            "AMR" : ["CLM", "MXL", "PEL", "PUR"],
            "EAS" : ["CDX", "CHB", "CHS", "JPT", "KHV"],
            "EUR" : ["CEU", "FIN", "GBR", "IBS", "TSI"],
            "SAS" : ["BEB", "GIH", "ITU", "PJL", "STU"]
    }

    # small loop to generate each line:
    def get_allele_frequencies(pop, alleles, data, pop_class):
        population_names = config.population_names

        try:
            popname = population_names[pop]
        except:
            popname = "NA"
            print pop
        # returning all frequencies:
        string = "\t<tr class=\"%s\"><td class=\"label\">%s (%s)</td>" %  (pop_class, popname, pop)
        for allele in alleles:

            freq = float()
            try:
                freq =  round(data[pop][allele], 4)

            except:
                freq = 0.00
            string += "<td>%s:%s</td>" %(allele, str(freq))

        string += "</tr>\n"
        return string;

    # Defining the header of the table:
    table += "<table class=\"populations\">\n\t<tr class=\"pop_header\"><td>Population (code)</td><td >Allele 1</td><td>Allele 2</td></tr>\n"
    table += get_allele_frequencies("ALL", alleles, pops, "pop_ALL")
    for big in population_groups.keys():
        table += get_allele_frequencies(big, alleles, pops, "pop_"+big)

        for small in population_groups[big]:
            table += get_allele_frequencies(small, alleles, pops, "pop_"+small)

    table += "</table>\n"
    return table

#
def draw_gene_table(all_genes, filename):
    '''
    This function creates a html table for the report.
    Input:
        all_genes: the output of the get_gene_list function.
        filename: the filename of the generated html file.

    Output:
        html formatted table string.
    '''
    html_string = "<table class=\"gene_list\">\n"
    html_string += "\t<tr style=\"font-weight: bold\"><td>Gene name</td><td>Ensembl ID</td><td>Biotype</td><td>Distance</td><td>Orientation</td></tr>\n"

    for distance in sorted(all_genes.keys()):

        for gene in all_genes[distance]:

            gene_name = gene["name"]

            try:
                description = gene["description"].split("[")[0]
            except:
                description = "-"
            gene_ID = gene["ID"]
            biotype = gene["biotype"]
            distance = gene["distance"]
            orientation = gene["orientation"]

            # formatting distance string:
            if distance > 1000:
                distance = distance / 1000
                distance = str(distance)+"kbp"
            else:
                distance = str(distance) + "bp"

            # selecting class:
            div_class = "other_gene"
            if orientation == "overlapping":
                div_class = "overlap"
            elif biotype == 'protein_coding':
                div_class = "protein_coding"

            # Generating the line:
            html_string += "\t<tr class=\"%s\">" % (div_class)
            a_tag = "<a href=\"./genes/%s.html\">%s</a>" %(gene_ID, gene_name)

            html_string += "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" %(a_tag, description, biotype, distance, orientation)


    html_string += "</table>"
    return html_string

#
def draw_consequence_table(VEP_data):
    '''
    This function draws a html table with the consequences.
    It handles only consequences with overlapping transcripts.
    '''

    if type(VEP_data) == str:
        html = '<div class="missing_data">%s</div>' % GTEX_data
        return html

    # Creating table header:
    table_row = "\t<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n"
    Ensembl_link = '<a href="http://grch37.ensembl.org/Homo_sapiens/Regulation/Summary?fdb=funcgen;rf=%s">%s</a>'
    APPRIS_link = '<a href=\"http://appris.bioinfo.cnio.es/#/database/id/homo_sapiens/%s?db=hg19">%s</a>'

    table = ""

    if type(VEP_data['transcript']) == str:
        table += '<div class="missing_data">%s</div>' % VEP_data['transcript']
    else:
        table += "The following transcripts overlap with the queried variation:<br>\n"

        table += "<table class=\"consequence_table\">\n"
        table += "\t<tr><td class=\"table_header\">Gene</td><td class=\"table_header\">Transcript ID</td><td class=\"table_header\">Impact</td><td class=\"table_header\">Consequence</td><td class=\"table_header\">Principial Isoform</td></tr>\n"

        # Filling data:
        for consequence in VEP_data["transcript"]:
            transcript_ID = Ensembl_link % (consequence['transcript_id'], consequence['transcript_id'])
            gene_ID = Ensembl_link % (consequence['gene_id'], consequence['gene_symbol'])
            APPRIS_field = APPRIS_link % (consequence['gene_id'], consequence['principal'])
            table +=  table_row % (gene_ID, transcript_ID , consequence['impact'], ",<br>".join(consequence['consequence']), APPRIS_field)
        table += "</table><br>\n"

    if type(VEP_data['regulatory']) == str:
        table += '<div class="missing_data">%s</div>' % VEP_data['regulatory']
    else:
        table += "The following regulatory features overlap with the queried variation:<br>\n"
        table += "<table class=\"consequence_table\">\n"
        table += "\t<tr><td class=\"table_header\">Biotype</td><td class=\"table_header\">Regulatory feature ID</td><td class=\"table_header\">Impact</td><td class=\"table_header\">Consequence</td></tr>\n"
        table_row = "\t<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n"
        # Filling data:
        for consequence in VEP_data["regulatory"]:

            # If any of the values are missing, we have to add fake elements:
            try:
                biotype = consequence['biotype']
            except:
                biotype = "NA"

            try:
                impact = consequence['impact']
            except:
                impact = "NA"

            try:
                cons = ",<br>".join(consequence['consequence'])
            except:
                cons = "NA"

            ID = Ensembl_link % (consequence['regulatory_ID'], consequence['regulatory_ID'])
            table +=  table_row % (biotype, ID, impact, cons)
        table += "</table><br>"

    return table

#
def draw_pubmed_table(pubmed_data):
    '''
    This function creates a formatted table with all the pubmed entries where the rsID was mentioned in anywhere in the text.
    '''

    table = "<table class=\"consequence_table\">\n"

    for ID in pubmed_data.values():
        table += "\t<tr><td><b>%s</b> <i>et al.</i> %s (%s): <a href=\"%s\">%s</a><br></td></tr>\n" %(ID['firstAuthor'],ID['journal'],ID['year'],ID['URL'],ID['title'],)

    table += "</table>"

    return table.encode('utf-8')

# This function was developed to create html formatted tables from the returned regulatory data:
def draw_regulation(regulatory_data):
    '''
    This function creates html formatted tables of the returned regulatory features.

    input: output of the regulatory function.
    output: html formatted string.
    '''

    header = {
        'Transcription' : "<br>Transcription factor binding sites:\n",
        'Open' : "<br>Open chromatin:\n",
        'Polymerase' : "<br>Polymerase binding:\n",
        'Histone' : "<br>Histone modifications:\n"
    }

    html = "<br>"

    for feature in regulatory_data.keys():
        html += header[feature]

        # now generating the html formatted table:
        html += "<table id=\"regulation\" class=\"regulation\">\n"
        for reg_type in regulatory_data[feature].keys():
            html += "\t<tr><td style=\"width:75px;font-weight:bold\">%s</td><td>%s</td></tr>\n" %(
                    reg_type, "/".join(regulatory_data[feature][reg_type]))

        html += "</table>"

    return (html)

# This function draws a html formatted table of the ExAC frequency data:
def draw_ExAC_table(Exac_parsed, ref):
    '''
    This function creates a html formatted table with the ExAC data previously generated

    Input data: dictionary generated by get_ExAC

    '''
    table ='<br>\n<div id="ExAC_data" class="general" >Frequencies published by the Exome Aggregation Consortium:</div>\n'

    # If the returned value is a string we report the returned value:
    if CheckReturnedValue(Exac_parsed): return CheckReturnedValue(Exac_parsed)

    ExAC_pops = ['NFE','FIN','AFR','AMR','EAS','SAS','OTH','ALL']
    ExAC_pops_explained ={
        'AFR' : "African",
        'ALL' : "Total",
        'FIN' : "European (Finnish)",
        'AMR' : "Latino",
        'EAS' : "East Asian",
        'SAS' : "Sout Asian",
        'OTH' : "Other",
        'NFE' : "European (Non-Finnish)"
    }

    # Writing header:
    table += "<table class=\"ExAC_pop\">\n"
    table += "\t<tr class=\"pop_header\"><td>Population</td>"
    for allele in Exac_parsed.keys():
        table += "<td>Alt: %s, freq. (count)</td>" % allele
    table += "<td>Ref: %s freq.</td></tr>\n" %ref

    # Looping through all populations:
    for pop in ExAC_pops:
        #print ExAC_pops_explained[pop]
        table += "\t<tr><td class=\"table_header\">%s</td>" % ExAC_pops_explained[pop]
        ref_freq = 1
        for allele in Exac_parsed.keys():
            count = Exac_parsed[allele][pop]["count"]
            freq = round(float(Exac_parsed[allele][pop]["frequency"]),5)
            ref_freq -= float(Exac_parsed[allele][pop]["frequency"])
            table += "<td>%s (%s)</td>" % (freq, count)

        # Adding reference freq:
        table += "<td>%s</td></tr>\n" % (round(ref_freq, 5))

    table += "</table>\n"
    return table

# Creates a html formatted table of uk10k allele frequencies:
def draw_UK10K_table(ukData):
    '''
    Input: dictionary create by get_UK10K_frequencies()
    Input structure:
    {'Allele_counts': {
        'A': 7446,
        'G': 116},
     'Allele_frequencies': {
        'A': 0.98466014281936,
        'G': 0.015339857180640043},
     'Genotype_counts': {
        'AA': '3665 (3665.89)',
        'GA': '116 (114.22)',
        'GG': '0 (0.89)'}}

    If there was some problem, a sting is returned with the status.

    Output: html formatted table.
    '''

    UK10K_table ='<br>\n<div class="general">Allele frequencies and allele counts in the UK10K data:</div>\n'

    # If the returned value is a string we report and return with a smaller output:
    if type(ukData) is str:
        UK10K_table += '<div class="missing_data">%s</div>\n' %(ukData)
        return UK10K_table
    else:
        pass

    # Initializing table:
    UK10K_table += '<table class="UK10k_table">\n'
    # Fill table with data:
    for allele in ukData['Allele_frequencies']:
        allele_freq = ukData['Allele_frequencies'][allele]
        allele_count = ukData['Allele_counts'][allele]

        UK10K_table += '\t<tr><td>"%s" allele frequency (count)</td><td>%s (%s)</td></tr>\n' %(allele, allele_freq, allele_count)

    for genotype in ukData['Genotype_counts']:
        genotype_count = ukData['Genotype_counts'][genotype]

        UK10K_table += '\t<tr><td>"%s" genotype count (expected)</td><td>%s</td></tr>\n' %(genotype, genotype_count)

    # Closing table:
    UK10K_table += "</table>\n"
    # Adding comment:
    UK10K_table += '<div class="table_comment">\n\t%s Expected genotype counts were calculated assuming HW equlibrium.\n</div>\n' % (config.UK10K_VERSION)

    return UK10K_table

# Drawing tables of genes affected by the queried variation:
def draw_GTEx_eQTL_table(GTEX_data):
    '''
    This function creates an html formatted string based on the previously retrieved GTEx data.
    Input: output of the get_GTEx_genes function that can contain error messages as well.
    Output: html formatted string that can be used directly to generate the variation html file.
    '''

    # Check if the returned value is a dictionary or a string:
    if CheckReturnedValue(GTEX_data): return CheckReturnedValue(GTEX_data)

    # If a dictionary, initialize table:
    table = '<table class="GTEX_eQTL">\n'
    table += '\t<tr><td>Tissue</td><td>p-value</td><td>beta (sd)</td><td>Gene name</td><td>Distance from TSS</td></tr>\n'

    # Get basic information:
    Tissue_count = len(GTEX_data.keys())
    Tissues = GTEX_data.keys()
    Gene_count = 0
    Genes = []

    for tissue in GTEX_data.values():
        for association in tissue:
            #var = variation % (association["rsID"],association["rsID"])
            gen = config.ENSEMBL_GENE % (association["gene_ID"], association["gene_name"])
            table += "\t<tr><td>%s</td><td>%s</td><td>%s (%s)</td><td>%s</td><td>%s</td></tr>\n" % (
                association["tissue"], association["pvalue"], association["beta"],association["sd"],  gen, association["distance"])

            # Adding gene to gene collection:
            if not association["gene_name"] in Genes:
                Genes.append(association["gene_name"])
                Gene_count += 1

    table += '</table>\n'

    # Generating a summary of the returned data:
    Report = "The expression of "
    if Gene_count > 1: Report += "%s genes (<b>%s</b>) were found to be affected by %s" %(Gene_count, ", ".join(Genes), association["rsID"])
    else : Report += "%s gene was found to be affected by %s" %(Genes[0], association["rsID"])
    if Tissue_count > 1: Report += " in %s tissues " %( Tissue_count)
    Report += '(' + config.GTEX_LINK % association["rsID"] + '):'

    # Adding report and table together in a html formatted string:
    html = '<div class="GTEx variation" onclick="toggle()" style="margin-left:10px; margin-bottom:2px">%s\n' % Report
    html += '<div class="GXA_hidden" style="display: none;">\n%s\n</div></div>\n' %table

    return html

# Processes tables with phenotype data:
def draw_phenotpye_table_var(df):
    '''
    Ths function creates an html formatted table based on the submitted dataframe
    input: dataframe
    output: html formatted table

    If string is submitted we just print the text.
    '''
    ## Check if the returned value is a string or a dataframe:
    if type(df) == str or type(df) == unicode:
        html = '<div class="missing_data">%s</div>' % df
        return html

    def _process_phenotype_row(row):
        '''
        Helper function for generation tables
        '''
        URL_var = 'http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=%s' % row["rsID"]
        if row["distance"] > 1000:
            distance = row["distance"] / 1000
            distance = str(distance)+"kbp"
        else:
            distance = str(row["distance"])+"bp"
        return '\t<tr><td><a href="%s">%s</a></td>' \
                '<td>%s</td><td>%s</td><td><a href="%s">%s</a></td></tr>'  % (
                        URL_var, row["rsID"], row["consequence"], distance,
                        row['URL'], row["phenotype"]
                    )

    # Now we have to generate table header:
    html = '<table class="phenotypes">\n\t<tr class="pop_header"><td style="width:10%">rsID</td><td>Consequence</td><td>Distance</td>'\
            '<td>Phenotype</td></tr>\n'
    html += "\n".join(df.sort_values(by="distance").apply(_process_phenotype_row, axis=1))
    html += "\n</tbody></table>\n\n"
    return html