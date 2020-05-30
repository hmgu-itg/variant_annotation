'''
This package is a collection of functions shared between variations and genes

Version: 0.1 Last modified: 2016.02.24

Input and output of functions are explained at each function
'''

import requests
import datetime
import sys
import os
import re
import jinja2
import json 
import pandas as pd


# Footer generator:
def FooterGenerator(version):
    '''
    This function generates a html formatted footer that can be inserted into
    the htmls

    Input: version
    Output: simple string
    '''

    # Generating footer string:
    date_string = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    return "Generated by VarAnnoTool v.%s on: %s" %(version, date_string)

# Reading config file:
def ReadParameters(parameterFile):
    '''
    This function reads the parameter file, builds a hash from it, where the
    keys are the parameter names through wich other functions will be able to
    access the value of the parameters.

    Input: filename
    Oputput: dictionary
    '''

    # Checking if the parameterfile exists or not:
    if not os.path.isfile(parameterFile):
        sys.exit("Parameter file was not found: %s" % (parameterFile))

    commentChar = "#" # Lines in the parameter files starting with this character will be ignored
    separatorChar = "=" # This character will separate the keys and values in the parameter file
    ParameterDict = {} # Output

    # Opening the file:
    with open(parameterFile) as f:
        for line in f:
            line = line.strip()

            # Omitting empty lines:
            if not line:
                continue

            # Omitting lines which starts with comment character
            if line[0] == commentChar:
                continue

            # Parsing keys and values:
            m = re.match("(.+?)=[\"\'](.+)[\"\']", line)
            (key, value) = m.groups()
            ParameterDict[key.strip()] = value.strip()

    return ParameterDict

def restQuery(URL,qtype="get",timeout=None):
    func=None

    if qtype=="get":
        func=requests.get
    elif qtype=="post":
        func=requests.post
    else:
        print(str(datetime.datetime.now())+" : getQuery: query type ("+qtype+") has to be either \"get\" or \"post\"",file=sys.stderr)
        sys.stderr.flush()
        return None

    try:
        r = func(URL,headers={"Content-Type" : "application/json", "Accept" : "application/json"},timeout=timeout)
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


# REST submission:
def submit_REST (URL):
    '''
    This function submits the query to the REST server, and processes the output
    The returned a dictionary. If there is a problem, the script prints out the failed URL.

    Iniput: valid URL for REST query
    Output:

    '''

    r = requests.get(URL, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        print "Failed query: ",URL
        r.raise_for_status()
        sys.exit()

    return r.json()

# Checking returned values:
def CheckReturnedValue(var):
    '''
    A check which is performed by each drawer function
    to know if the returned value is a dictionary or a string.

    Input: any kind of variable
    Output: 0 if the variable is not a string
    html formatted string if the submitted variable is a string -> ready to plot.


    usage:
    if CheckReturnedValue(var): return CheckReturnedValue(var)

    '''

    if type(var) == str or type(var) == unicode:
        return '<div class="missing_data">%s</div>' % (var)

# Html generator:
def draw_html(templateFile, data):
    '''
    This is the html templating script. It needs the template html file, and the printed data set (in dictionary format).

    Output is the html string that can be directly saved as html.
    '''

    # search path for template files:
    templateLoader = jinja2.FileSystemLoader( searchpath="/" )

    # An environment provides the data necessary to read and
    templateEnv = jinja2.Environment( loader=templateLoader )

    # Read the template file using the environment object.
    # This also constructs our Template object.
    template = templateEnv.get_template( templateFile )

    # Finally, process the template to produce our final text.
    outputText = template.render( data )

    return outputText

# Shared function to get variants with phenotype annotations:
def get_phenotype_vars(input_dict):
    '''
    This function will be used by both the var and the gene trenches

    input dictionary:
        {"chr" : <str>, "start" : int, "end" : int, "var" : int}
    if var is given the output will contain the distance from that position.

    output: pandas dataframe with the following columns:
    'SNPID', 'consequence', 'distance', 'genes', 'phenotype', 'rsID', 'source'
    (Distance is in bp, and is an absolute value!)
    '''

    # Do we have the position?
    pos = input_dict["var"] if "var" in input_dict else 0

    # Let's check if the start is below zero?
    if input_dict["start"] < 0 : input_dict["start"] = 0

    # Now we return all variants with phenotypes:
    URL = "http://grch37.rest.ensembl.org/overlap/region/human/"\
            "%s:%d-%d?feature=variation;variant_set=ph_variants;content-type=application/json" %(
                input_dict["chr"], input_dict["start"], input_dict["end"])

    # Return all variants with phenotypes and overlapping with the given region.
    variants = submit_REST(URL)

    # Let's test if we have something:
    if len(variants) == 0: return "[Warning] No variants with phenotypes were found in this region."

    # Now extract all rsIDs from the dictionary:
    rsIDs = []
    for var in variants:
        rsIDs.append(var["id"])

    # Now formulate REST post
    variantPOSTURL = "http://grch37.rest.ensembl.org/variation/homo_sapiens?phenotypes=1"
    r = requests.post(variantPOSTURL,
                      headers={ "Content-Type" : "application/json",
                                "Accept" : "application/json"},
                      data=json.dumps({"ids" : rsIDs}))
    # Check output:
    if not r.ok:
        return "[Warning] Phenotypes could not be returned!"

    # process returned data:
    decoded = r.json()
    varsWithPhenotyes = []
    for rsID in decoded:
        for phenotype in decoded[rsID]["phenotypes"]:

            # Generate phenotype link based on source:
            if phenotype["source"] == "ClinVar":
                link = "https://www.ncbi.nlm.nih.gov/clinvar/?term="+rsID
            elif phenotype["source"] == "HGMD-PUBLIC":
                link = "http://www.hgmd.cf.ac.uk/ac/gene.php?gene="+phenotype["genes"]
            elif "NHGRI-EBI" in phenotype["source"]:
                link = "https://www.ebi.ac.uk/gwas/search?query="+rsID
            elif phenotype["source"] == "OMIM":
                link = "https://www.omim.org/entry/"+phenotype["study"][4:]
            else:
                link = "http://grch37.ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;v=CM106680;vdb=variation"

            # Formulate return value:
            varsWithPhenotyes.append({"rsID": rsID,
                    "consequence" : decoded[rsID]["most_severe_consequence"],
                    "SNPID" : "chr%s:%s" %(
                        decoded[rsID]["mappings"][0]["seq_region_name"],
                        decoded[rsID]["mappings"][0]["start"]),
                    "phenotype" : phenotype["trait"],
                    "URL" : link,
                    "distance" : abs(pos - decoded[rsID]["mappings"][0]["start"])})

    # Create a dataframe:
    df = pd.DataFrame(varsWithPhenotyes)
    df.consequence = df.consequence.str.title().str.replace("_", " ")

    return df
