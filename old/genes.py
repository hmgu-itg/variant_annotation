'''
These functions were collected to annotate individual genes. Queries will retireve
information from various sources, so there will be various functions accessing
information through internet.

Package includes:
* submit_REST
* get_gene_info

version: v.2.1.0 Last modified: 2015.11.23

'''

# Libraries for configuration:
import config

# Shared functions:
from shared import *

import sys, os, re, requests, ast

# Importing date library:
import datetime

# Library to execute shell commands:
import subprocess

# list of functions that will create all tables and stuff:
from genes_draw import *

# Pandas dataframe to deal with certain data structures:
import pandas as pd

# list of functions that will create all tables and stuff:
#from variations_draw import *

# Upon loading this module, automatically generate gene folder:
os.system('mkdir -p ./genes')

### Combining all steps togeteher to perform all annotation:
def Annotate_gene(gene_ID):
    '''
    This function performs all gene annotation and saves all data into an html file into a folder called gene.
    input: gene_ID
    output: nothing.
    '''

    # Version:
    version = config.VERSION

    # Generating footer string:
    footer = FooterGenerator(version)

    ### Retrieving information:
    # Retrieving general information:
    general_info = get_gene_info(gene_ID)

    # Retrieving cross references:
    cross_refs = get_gene_xrefs(gene_ID)

    # Retrieving uniprot data:
    Uniprot_data = get_UNIPROT_data(cross_refs)

    # Downloading information from OMIM database:
    OMIM_annotation = get_OMIM(cross_refs)

    # Generate GWAS list:
    GWAS_signals = get_gwas_hits(general_info)

    # Retrieve data from gene expression atlas:
    (GXA_levels, GXA_data) = get_GXA(general_info)

    # Retrieving GTEx data:
    GTEx_data = get_GTEx_variations(general_info)

    # Get mouse phenotypes:
    mouse_pheno = get_mouse_phenotype(gene_ID)

    ### Drawing retrieved data:
    gene_details = draw_gene_info(general_info)

    # Drawing GO annotations:
    if len(cross_refs["GO"]) > 1:
        GO_table = draw_GO_terms(cross_refs, gene_ID)
    else:
        GO_table = ""

    # Drawing GXA table:
    GXA_table = draw_GXA_table(GXA_levels, general_info)

    # Drawing GXA heatmap:
    GXA_heatmap = draw_GXA_heatmap(GXA_data)

    # Drawing GWAS signals:
    GWAS_table = draw_gwas_gene_table(GWAS_signals)

    # Drawing Uniprot table:
    Uniprot_table = draw_Uniprot_table(Uniprot_data)

    # Drawing GTEx table:
    GTEx_table = draw_GTEx_variations(GTEx_data)

    # Drawing OMIM table:
    OMIM_text = draw_OMIM_annotation(OMIM_annotation)

    # Generate html table with mouse phenotypes:
    mouse_pheno_table = draw_mouse_phenotypes(mouse_pheno)

    gene_dict = { "gene_table" : gene_details.decode('utf-8'),
                  "footer" : footer.decode('utf-8'),
                  "go_data": GO_table.decode('utf-8'),
                  "Uniprot_data" : Uniprot_table.decode('utf-8'),
                  "OMIM_annot" : OMIM_text.decode('utf-8'),
                  "GWAS" : GWAS_table.decode('utf-8'),
                  "GXA_table" : GXA_table.decode('utf-8'),
                  "GXA_heatmap" : GXA_heatmap.decode('utf-8'),
                  "GTExVariants" : GTEx_table,
                  "mouse_pheno" : mouse_pheno_table}

    # Gene specific template file:
    templateFile = config.GENE_TEMPLATE

    # Compiling html file:
    html = draw_html(templateFile, gene_dict)

    # Saving html data into a file. We need only the gene ID for the name:
    filename = "./genes/%s.html" % gene_ID
    f = open(filename, 'w')
    f.write(html.encode("utf8"))
    f.close()

# Download gene information from the Ensembl:
def get_gene_info (ID):
    '''
    This function retrieve every relevant information from Ensemble that we can
    access through the REST API

    INPUT: Ensembl stable ID
    OUTPUT: dictionary with retrieved information

    '''

    URL = config.REST_URL + "/lookup/id/%s" % ID

    response = submit_REST(URL)

    # Selecting info:
    gene_info = {}
    gene_info["source"] = response["source"]
    gene_info["id"] = response["id"]
    gene_info["start"] = response["start"]
    gene_info["end"] = response["end"]
    gene_info["assembly_name"] = response["assembly_name"]
    try:
        gene_info["description"] = response["description"].split("[")[0]
    except:
        gene_info["description"] = "NA"
    gene_info["name"] = response["display_name"]
    gene_info["type"] = response["biotype"]
    gene_info["strand"] = response["strand"]
    gene_info["chromosome"] = response["seq_region_name"]

    return gene_info

# Download cross references of a gene based on Ensembl annotation:
def get_gene_xrefs (ID):
    '''
    This function retrieve every relevant cross-references from Ensemble that we can
    access through the REST API
    '''

    URL = config.REST_URL + "/xrefs/id/%s?content-type=application/json;all_levels=1" % ID

    response = submit_REST(URL)

    # The following sources are implemented:
    xrefs = {
        "MIM disease" : [],
        "MIM gene" : [],
        "GO" : [],
        "GOSlim GOA" : [],
        "UniProtKB/Swiss-Prot" : [],
        "Human Protein Atlas" : [],
        "ChEMBL" : [],
    }

    # Now extracting cross references for all important sources:
    for xref in response:
        db =  xref['db_display_name']
        try:
            xrefs[db].append([xref["description"], xref["primary_id"]])
        except:
            continue

    return xrefs

# Retrieve data from the OMIM database through the API service:
def get_OMIM (cross_refs):
    '''
    Required information:
        - list of OMIM ID-s

    Retrieved information:
        - OMIM entry name
        - text field name
        - text field (that includes references as well)

    Retrieved data structure (one dict for all queried entries):
    dict{
        OMIM_ID :{
            title : OMIM_preferredTitle
            text : {
                OMIM_textSectionTitle : OMIM_textSectionContent
    }}}

    More information how the OMIM API works see: http://omim.org/help/api
    The API key will expire in 2016.11.11, at that time a new key have to be required.
    '''

    # Extracting OMIM ID from the cross refereces data:
    MIM_IDs = [x[1] for x in cross_refs["MIM disease"]]
    MIM_IDs.extend([x[1] for x in cross_refs["MIM gene"]])

    # The function returns at this point if there is not MIM IDs in the crossref
    if len(MIM_IDs) == 0:
        return "No OMIM entry for this gene."

    # Constructing query string:
    URL = config.OMIM_URL

    for ID in MIM_IDs:
        URL += 'mimNumber=%s&' % ID

    # Retrieving the following fields:
    URL += 'include=text&'
    URL += 'include=allelicVariantList&'
    URL += 'include=referenceList&'

    # Adding API key:
    URL += config.OMIM_APIKEY

    # Required format is pyton (although it is a string indeed)ormat=python
    URL += "format=python&"

    # Downloading page:
    page = requests.get(URL)

    # Reconstructingh dictionary from the returned string:
    OMIM_entry  = ast.literal_eval(page.content)

    # hash to fill in:
    OMIM_data = {}

    # Parsing hash:
    for entry in OMIM_entry['omim']['entryList']:
        # GEt OMIM ID:
        ID = entry["entry"]["mimNumber"]
        OMIM_data[ID] = {}

        # Get OMIM name:
        OMIM_data[ID]['title'] = entry["entry"]["titles"]['preferredTitle']

        # Get OMIM text:
        OMIM_data[ID]['text'] = {}
        for fields in entry['entry']["textSectionList"]:
            OMIM_data[ID]['text'][fields['textSection']['textSectionTitle']] = fields['textSection']['textSectionContent']

        # now we have to parse allelic variants:
        # print stuff['omim']['entryList'][0]['entry']['allelicVariantList'][0]['allelicVariant'].keys()
        # ['status', 'name', 'dbSnps', 'text', 'mutations', 'number', 'alternativeNames', 'clinvarAccessions']

        try:
            OMIM_data[ID]['variations'] = {}
            for variations in entry['entry']['allelicVariantList']:
                OMIM_data[ID]['variations'][variations['allelicVariant']['dbSnps']] = variations['allelicVariant']['text']
        except:
            continue

    return OMIM_data

# Finding already known GWAS associations for the gene:
def get_gwas_hits (general_info):
    '''
    This function performs a query for gwas signals in a local gwas datafile.
    All reported and mapped genes will give a hit.

    Keep in mind that this is a local datafile, so it has to be kept up to date to make sure
    it's still reliable
    '''

    # Extracting name from the general info:
    gene_name = general_info["name"]

    # path to gwas file:
    gwas_file = config.GWAS_FILE

    query = "grep -iw %s %s" % (gene_name, gwas_file)

    # gwas data is collected in this dictionary:
    gwas_data = []

    output = subprocess.Popen(query, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in output.stdout.readlines():
        fields = line.split("\t")

        # update dictionary for each hits:
        gwas_data.append({
            "rsID": fields[20],
            "p-value": fields[27],
            "Consequence": fields[24],
            "Reported genes": fields[13],
            "Mapped genes": fields[14],
            "trait": fields[7],
            "PMID": fields[1],
            "URL": fields[5],
            "Author": fields[2]
        })

    if not gwas_data:
        gwas_data = "No GWAS signal associated with the gene has been found."

    return gwas_data

# Retrieve uniprot data:
def get_UNIPROT_data(cross_refs):
    '''
    The cross references are taken as input. It might contain none, one or
    multiple references to the uniprot database. The function loops through all
    and returns a dictionary with all unique protein names as keys.

    Retrieved fields:
        Function, Subunit, Developmental stage,
        Tissue specificity, Catalytic activity,
        Disruption phenotype, Subcellular localiztion,
        Disease

    The list can be extended if need....
    '''

    def  __UNIPROT_DOWNLOAD__(ID):

        # Constructing URL
        URL = config.UNIPROT_URL
        URL += "?query=id:" + ID

        # Retrieved fields:
        URL += "&columns=id%2Ccomment%28FUNCTION%29%2C" # Function
        URL += "comment%28SUBUNIT%29%2C" # Subunit
        URL += "comment%28DEVELOPMENTAL%20STAGE%29%2C" # Developmental stage
        URL += "comment%28TISSUE%20SPECIFICITY%29%2C" # Tissue specificity
        URL += "comment%28CATALYTIC%20ACTIVITY%29%2C" # Catalytic activity
        URL += "comment%28DISRUPTION%20PHENOTYPE%29%2C" # Disruption phenotype
        URL += "comment%28SUBCELLULAR%20LOCATION%29%2C" # Subcellular localization
        URL += "comment%28DISEASE%29%2C" # Disease

        # Format:
        URL += "entry%20name&format=tab" # tab delimited format returned

        # Retrieve data from swissprot:
        r = requests.get(URL)

        # Processing data:
        fields = r.content.split("\n")[1].split("\t")

        # Filling data:
        UniprotData = {
            "Name" : fields[9],
            "Disease" : fields[8],
            "Function": fields[1],
            "Entry" : fields[0],
            "Subunit" : fields[2],
            "Phenotype" : fields[6],
            "localization" : fields[7],
            "Tissue" : fields[4],
            "Development" : fields[3]
        }

        return UniprotData

    # Looping through all cross references pointing to Uniprot:
    Independent_IDs = {}
    for uniprotID in cross_refs["UniProtKB/Swiss-Prot"]:
        Independent_IDs[uniprotID[0]] = uniprotID[1]

    # Annotate all unique proteins:
    annotated_protein  = {}
    for ID in Independent_IDs.keys():
        annotated_protein[ID] = __UNIPROT_DOWNLOAD__(Independent_IDs[ID])

    # Returning annotation:
    return (annotated_protein)

# Retrievig data from the EBI's gene expression atlas:
def get_GXA(general_info):
    '''
    This function downloads information from the gene expression atlas of EBI.
    Input: general info data structure, from which the Ensembl stable ID will be used
    Output: pandas core series with the 10 most highly expressing tissues.
    '''
    # Extracting name from the general info:
    gene_ID = general_info["id"]

    # the URL pointing to the gene expression atlas:
    URL = config.GXA_URL % gene_ID

    # The downloaded data will be read as a pandas dataframe:
    try:
        df = pd.read_csv(URL, sep='\t', comment="#")
    except:
        return ("No gene expression information!", "No gene information.")

    cleanExperiments = []
    for exp in df.Experiment.tolist():
        m = re.search('Tissues -\s+\d*(.+)', exp)
        cleanExperiments.append(m.groups()[0].strip())

    df.Experiment = cleanExperiments

    IndexNames = df.Experiment.tolist()

    # These are the main sources that we focus on:
    sources = ["GTEx", "FANTOM5 project - adult", "FANTOM5 project - fetal"]
    Levels = {}

    # Extracting top 10 tissues for each sources:
    for source in sources:
        try:
            fieldName = filter(lambda x : source in x, IndexNames)[0]
            Index = IndexNames.index(fieldName)

            # Excising the GTEx data from dataframe:
            GTEx = df.ix[Index]

            GTEx = GTEx.sort_values(ascending=False, na_position='last')[0:11]
            Levels[fieldName] = []
            for tissue, value in GTEx.iteritems():
                try:
                    if float(value) > 0 : Levels[fieldName].append([tissue, value])
                except:
                    pass
        except:
            pass

    return (Levels, df)

# checking which variations effect the expression of a given gene:
def get_GTEx_variations(general_info):
    '''
    This function returns a dictionary with variations that have an effect on the
    expression of the gene. Each key of the dictionary is a variation, and the
    values are dictionaries that includes the observed tissues, p-values.

    Always check for GTEx website! -> A link should be include for the gene.
    '''

    # In the upcoming version this has to be read from the configuration file:
    GTEx_file = config.GTEX_FILE_GENES

    # Check if datafile is exists, and return 0 if not:
    if not os.path.isfile(GTEx_file):
        return "GTEx datafile was not found in this location: %s" % GTEx_file

    # Extract data for bedtools query:
    chromosome = general_info["chromosome"]
    start = general_info["start"]
    end = general_info["end"]
    ID = general_info["id"]

    # submit bcftools query:
    query = "bash -O extglob -c \'/software/hgi/pkglocal/tabix-git-1ae158a/bin/tabix %s chr%s:%s-%s | grep %s\'" %(GTEx_file, chromosome, start, end, ID)
    output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    GTEx_data = dict()
    for line in output.stdout.readlines():
        fields = line.strip().split("\t")
        rsID = fields[4]
        tissue = fields[3]
        pvalue = fields[11]
        distance = fields[8]
        gene_name = fields[5]

        # Initialize hit:
        if not rsID in GTEx_data:
            GTEx_data[rsID] = []

        # Add fields to data:
        GTEx_data[rsID].append({
                "rsID" : rsID,
                "tissue" : tissue,
                "distance" : distance,
                "pvalue" : pvalue,
                "gene_ID": general_info["id"],
                "gene_name" : gene_name
            })
    if len (GTEx_data) == 0:
        GTEx_data = "No eQTL signal was found for gene %s" %(general_info["name"])

    return GTEx_data

def _get_mouse_ID (human_ID):
    '''Looking up mouse gene ID of a given human gene ID'''

    URL = "http://grch37.rest.ensembl.org/homology/id/%s?content-type=application/json&target_species=mouse&aligned=0&sequence=none&type=orthologues" % human_ID
    data = submit_REST(URL)

    # If there is a problem:
    if "error" in data: return data["error"]

    # Now from the returned data we have to pick all possible homolgue IDs:
    mouse_IDs = {}
    if len(data["data"]) == 0: return "[Info] No mouse cross-ref for this gene!"
    if len(data["data"][0]["homologies"] ) == 0: return "[Info] No mouse cross-ref for this gene!"

    for homolog in data["data"][0]["homologies"]:
        try:
            URL = "http://grch37.rest.ensembl.org/lookup/id/%s?content-type=application/json;expand=0" % homolog["target"]["id"]
            data = submit_REST(URL)
            mouse_IDs[data["id"]] = [data["description"].split(" [")[0],
                              data["display_name"]]
        except:
            continue
    if len(mouse_IDs) == 0: mouse_IDs = [("-","-","-")]
    return mouse_IDs

def _get_MGI_ID (mouse_ID):
    '''Looking up MGI cross reference for a given mouse gene ID.'''

    URL = "http://grch37.rest.ensembl.org/xrefs/id/%s?content-type=application/json&external_db=MGI" % mouse_ID
    data = submit_REST(URL)

    # Now from the returned data we have to pick all possible homolgue IDs:
    for d in data:
        try:
            return d["primary_id"]
        except:
            continue
    return "[Info] MGI ID was not found"

def _get_MGI_phenotypes (MGI_ID):
    '''returning phenotype information stored on http://www.informatics.jax.org/ '''

    # Download whole site:
    URL = "http://www.informatics.jax.org/allele/report.txt?markerId=%s" % MGI_ID

    # Return all associated allele information:
    try:
        df = pd.read_csv(URL, sep="\t")

        # Drop unnecessary columns:
        df.drop([ u'Allele Symbol', u'Chromosome',  u'Synonyms',u'Allele Attributes', u'Transmission', u'Unnamed: 10'], axis=1, inplace=True)
        df.columns = [u'Allele_ID', u'Allele_name', u'Allele_type',
           u'Phenotypes',
           u'Human_disease']
        return df
    except:
        return "[Info] no phenotype was found!"

def get_mouse_phenotype (gene_id):
    '''returning mouse phenotype given human gene ID'''
    # Returning all mouse homologue IDs:
    mouse_gene_IDs = _get_mouse_ID(gene_id)

    # Checking if the returned data is a string, don't bother with the rest:
    if isinstance(mouse_gene_IDs, basestring): return mouse_gene_IDs

    # Looping through all mouse gene IDs, and get the MGI ID:
    MGI_IDs = {}
    for mouse_gene_ID in mouse_gene_IDs:
        MGI_IDs[mouse_gene_ID] = _get_MGI_ID(mouse_gene_ID)

    # Once we have all the MGI identifiers, we retrieve all the phenotypes:
    full_dataframe = ""
    for mouse_id, mgi_id in MGI_IDs.iteritems():
        df = _get_MGI_phenotypes(mgi_id)

        # Adding extra columns for the record:
        df["mouse_gene_ID"] = mouse_id
        df["MGI_ID"] = mgi_id
        df["mouse_gene_name"] = mouse_gene_IDs[mouse_id][1]
        df["mouse_gene_description"] = mouse_gene_IDs[mouse_id][0]
        if isinstance(full_dataframe, basestring):
            full_dataframe = df
        else:
            full_dataframe = pd.merge(full_dataframe, df, how="outer")

    return full_dataframe

