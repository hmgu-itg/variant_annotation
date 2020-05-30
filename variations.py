'''
These functions were collected to annotate individual variations.

version: v.2.1.0 Last modified: 2015.11.23

Package includes:

'''

# Libraries for configuration:
import config

# Shared functions:
from shared import *

# library for simple web queries:
import requests, sys, os

# libray for parsing htlm files:
from lxml import html

# library for generate python data structure from string:
import ast

# Library to execute shell commands:
import subprocess

def get_GWAVA_score(variant_data):
    '''
    This function calculates the gerp and gwava scores for a given variation.
    Of course it calls the gwava on farm, and parses the result.
    '''

    # Using variant_data for input:
    bed_string = "chr%s\t%s\t%s\t%s\n" % (variant_data["chromosome"], int(variant_data["start"]) - 1, variant_data["start"],variant_data["rsID"])

    # temporary input filename:
    filename = "/tmp/chr%s.bed" % variant_data["location"].replace(":", "_")

    # temporary output filename with gwava annotation:
    annot_filename = filename+"_ann"

    # temporary output filename with gwava prediction:
    gwava_filename = filename+"_gwava"

    # Saving primary input file:
    f = open( filename, 'w')
    f.write(bed_string)
    f.close()

    # now we have to run gwava:
    GWAVA_dir = config.GWAVA_DIR

    #query = "bash -O extglob -c \'python %s/src/gwava_annotate.py %s %s\'" %(GWAVA_dir, filename, annot_filename)
    query = "python %s/src/gwava_annotate.py %s %s" %(GWAVA_dir, filename, annot_filename)

    PATH = "/software/hgi/pkglocal/samtools-1.2/bin:/software/hgi/pkglocal/vcftools-0.1.11/bin:/software/hgi/pkglocal/tabix-git-1ae158a/bin:/software/hgi/pkglocal/bcftools-1.2/bin:/nfs/team144/software/ensembl-releases/75/ensembl-tools/scripts/variant_effect_predictor:/nfs/team144/software/bedtools2/bin:/nfs/team144/software/scripts:/nfs/users/nfs_d/ds26/bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/software/bin"

    # Submit query:
    output = subprocess.Popen(query.strip(),
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              env={'GWAVA_DIR': GWAVA_dir,
                                   "PYTHONPATH": "/nfs/team144/software/anaconda/lib/python2.7/site-packages",
                                   "PATH": PATH}).wait()

    # Once the annotation run is completed, let's run the prediction:
    query = "python %s/src/gwava.py tss %s %s" %(GWAVA_dir, annot_filename, gwava_filename)

    # Submit query:
    output = subprocess.Popen(query.strip(),
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              env={'GWAVA_DIR': GWAVA_dir,
                                   'PYTHONPATH': '/nfs/team144/software/anaconda/lib/python2.7/site-packages',
                                   'PATH': PATH}).wait()

    # Once the queries are returned, we have to parse the output:
    # From the annotation file, we retrive 147 - avg_gerp, 148 - gerp

    # Reading annotation file:
    for line in open(annot_filename, 'r'):
        if "start" in line: continue

        avg_gerp = line.split(",")[147]
        gerp = line.split(",")[148]

    # Reading prediction file:
    for line in open(gwava_filename, 'r'):
        gwava = line.strip().split("\t")[4]

    # Updating variant
    variant_data["avg_gerp"] = round(float(avg_gerp), 4)
    variant_data["gerp"] = gerp
    variant_data["gwava"] = gwava

    # removing temporary files:
    os.remove(annot_filename)
    os.remove(filename)
    os.remove(gwava_filename)

    return variant_data

# This function retrieves a list of gwas signals around the variation.
def get_gwas_hits_position(variant_data):
    '''
    This function retrives a list of all gwas signals within 500kbp distance around a given position.

    input: chromosome, position

    output: [ {"rsID","SNPID","trait","p-value","PMID","distance"}]
    '''
    # Gwas bedfile to use:
    gwas_file = config.GWAS_FILE_VAR

    # Check if datafile is exists, and return 0 if not:
    if not os.path.isfile(gwas_file):
        return "[Warning] Local GWAS datafile was not found in this location: %s\n" % gwas_file

    # Generating extended regions around our variation for bedtool search:
    start = variant_data["start"] - 500000
    if start < 0: start = 0
    end = variant_data["end"] + 500000
    chromosome = variant_data["chromosome"]

    # Assembling intersectBed query:
    # This is a bit tricky, as the standard subprocess calling is sh not bash. So we have to specify we are using bash this time.
    query = "bash -O extglob -c \'/nfs/team144/software/bedtools2/bin/intersectBed -b <(echo -e \"chr%s\t%s\t%s\") -a %s\'" % (chromosome, start, end, gwas_file)

    # Call intersectBed:
    output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Process output:
    signals = []
    for line in output.stdout.readlines():
        # chr1	3669356	3669356	rs76597070	Obesity-related traits	CCDC27	CCDC27	23251661	3E-6	0.00259585
        fields = line.strip().split("\t")

        # Calculating distance from our variation:
        distance = int(fields[1]) - int(variant_data["start"])

        if abs(distance) > 1000:
            distance = str(distance / 1000)+"kbp"
        else:
            distance = str(distance) + "bp"

        signals.append({
            "rsID"  : fields[3],
            "SNPID" : fields[0]+":"+fields[1],
            "trait" : fields[4],
            "p-value"  : fields[8],
            "PMID"     : fields[7],
            "distance" : distance
        })

    # return string if nothing has been found:
    if len(signals) == 0:
        signals = "No GWAS signals has been found around the variation."

    return signals

# this function is designed to return a chunk of the reference genome
def get_refSeq(chromosome, start, end):
    '''
    this function retrieves a certain chunk of the reference genome, based on the submitted
    chromosome, end and start positions.

    The returned value is the plain sequence.
    '''

    URL = config.REST_URL + "/sequence/region/human/%s:%s..%s:1" % (chromosome, start, end)
    response = submit_REST(URL)
    return response["seq"]

# Retiveing data from the variant effect predictor server.
def get_VEP_data(variant_data):
    '''
    This function returns the predicted effect based on chromosome, position and the alternate allele.

    Input: variant data given by the get_variant_info() function.

    Output: two dictionaries returned: list of variant effect, and an updated variant data
    '''

    # Parsing input data
    chromosome = variant_data["chromosome"]
    start = variant_data["start"]
    end = variant_data["end"]
    if end < start: start,end = end, start

    # OK, at first we check the reference allele, and we keep the allele which is not matches the reference.
    # There is no error checking here. If none of the alleles matches the reference allele, it won't show up.
    # Maybe it should be fixed later.
    refseq = get_refSeq(chromosome, start, end)
    allele = variant_data["allele_string"].split("/")[1]
    if allele == refseq: allele = variant_data["allele_string"].split("/")[0]

    URL = config.REST_URL + "/vep/human/region/%s:%s-%s:1/%s?content-type=application/json" % \
        (chromosome, start, end, allele)

    VEP = submit_REST(URL)

    # The most severe Sift and Polyphen scores will be retrieved:
    sift = []
    polyphen = []

    # Big dataset for all transcripts and regulatory features will be retrieved:
    VEP_data = {'transcript' : [],
                'regulatory' : []}

    # We retrieve the amino acid changes as well:
    aminoacid = ""
    codon = ""

    # If there is no consequence in the variant data then, we add: <-It does not work....
    variant_data["consequence"] = VEP[0]['most_severe_consequence']

    # Gathering transcript consequences:
    if 'transcript_consequences' in VEP[0].keys():
        VEP_data['transcript']=[]
        for transcript in VEP[0]['transcript_consequences']:
            VEP_data['transcript'].append({
                "impact" : transcript["impact"],
                "gene_symbol" : transcript["gene_symbol"],
                "gene_id": transcript["gene_id"],
                "transcript_id": transcript["transcript_id"],
                "consequence" : transcript["consequence_terms"]
            })

            # Testing polyphen score:
            if "polyphen_score" in transcript.keys():
                try:
                    if transcript["polyphen_score"] > polyphen[0]:
                        polyphen[0] = transcript["polyphen_score"]
                        polyphen[1] = transcript["polyphen_prediction"]
                except:
                        polyphen.append(transcript["polyphen_score"])
                        polyphen.append(transcript["polyphen_prediction"])

            # Testing sift score:
            if "sift_score" in transcript.keys():
                try:
                    if transcript["sift_score"] < sift[0]:
                        sift[0] = transcript["sift_score"]
                        sift[1] = transcript["sift_prediction"]
                except:
                        sift.append(transcript["sift_score"])
                        sift.append(transcript["sift_prediction"])

            # Testing amino acid substitution:
            if VEP[0]["most_severe_consequence"] in transcript["consequence_terms"] and "amino_acids" in transcript.keys():
                aminoacid = transcript["amino_acids"].replace("/", str(transcript["protein_start"]))
                codon = transcript["codons"]


    # Gathering regulatory consequences:
    if 'regulatory_feature_consequences' in VEP[0].keys():
        VEP_data['regulatory']=[]
        for regulatory in VEP[0]['regulatory_feature_consequences']:
            VEP_data['regulatory'].append({
                "impact" : regulatory["impact"],
                "biotype" : regulatory["biotype"],
                "gene_id": "-",
                "regulatory_ID": regulatory["regulatory_feature_id"],
                "consequence" : regulatory["consequence_terms"]
            })

    # Updating variation data:
    variant_data["sift"] = sift
    variant_data["polyphen"] = polyphen
    variant_data["codon"] = codon
    variant_data["aminoacid"] = aminoacid

    # If there were overlapping transcripts, we have to check if the
    # get a list of all overlapping genes:
    gene_list = []
    for transcripts in VEP_data['transcript']:
        if not transcripts["gene_id"] in gene_list:
            gene_list.append(transcripts["gene_id"])

    # Now looping trugh all genes and get list of principal transcripts:
    for gene_id in gene_list:
        appris_data = get_principal_transcript(gene_id) # downloading appris data

        # Now looping trough all consequences to get the appris annotation:
        for transcripts in VEP_data['transcript']:
            if transcripts["gene_id"] in gene_id:
                try:
                    transcripts["principal"] = appris_data[transcripts["transcript_id"]]["reliability"]
                except:
                    transcripts["principal"] = "NA"

    # Checking if there is anything in the VEP data:
    if len(VEP_data['transcript']) == 0:
        VEP_data['transcript'] = "The queried variation does not overlap any transcripts."
    if len(VEP_data['regulatory']) == 0:
        VEP_data['regulatory'] = "The queried variation does not overlap any regulatory features."

    # Returning data
    return (VEP_data, variant_data)


def get_principal_transcript(gene_ID):
    '''
    This function returns transcript information from the {appris} webserver based on Ensembl gene ID
    The output is a dictionary with detailed information of all annotated transcript.

    Link: http://appris.bioinfo.cnio.es/

    Known issue: appris gives error message if the queried gene was not found in the database. As
    most non-protein coding genes are missing from the database I don't want to inform the user
    about this error. Therefore, if there was some problem with the query, it will be lost.
    '''

    #
    URL = "http://apprisws.bioinfo.cnio.es:80/rest/exporter/id/homo_sapiens/%s?format=json&db=hg19" % gene_ID;

    r = requests.get(URL)
    if not r.ok:
        print "{appris} Query was not successful for gene: %s\nURL: %s" % (gene_ID,URL)
        return ""

    decoded = r.json()
    principal = {}
    for transcript in decoded:
        if not transcript["transcript_id"] in principal:
            principal[transcript["transcript_id"]] = {
                "reliability" : "NA",
                "transcript_ID" : transcript["transcript_id"],
                "name" : "NA",
                "length_aa" : "NA",
                "length_na" : "NA",
                "type" : ""
            }
        try: principal[transcript["transcript_id"]]["reliability"] = transcript["reliability"];
        except: pass
        try: principal[transcript["transcript_id"]]["transcript_ID"] = transcript["transcript_id"];
        except: pass
        try: principal[transcript["transcript_id"]]["name"] = transcript["transcript_name"];
        except: pass
        try: principal[transcript["transcript_id"]]["length_aa"] = transcript["length_aa"];
        except: pass
        try: principal[transcript["transcript_id"]]["length_na"] = transcript["length_na"];
        except: pass
        try: principal[transcript["transcript_id"]]["type"] += transcript["type"] + "/";
        except: pass

    return principal

# This function downloads variation data from ensembl:
def get_variant_info(rsID, ref):
    '''
    This function returns the most important details of a variation given it's rsID.

    It also returns the 1000 Genome Phase 3 allele frequencies as well.

    Returned: pop_freq, variation_data
    '''

    URL = config.REST_URL + "/variation/human/%s?content-type=application/json;pops=1;phenotypes=1" %rsID

    variationData = submit_REST(URL)

    # If the user only provides rsID, we don't know which is hte reference allele:
    if ref == "NA":
        ref = variationData["ancestral_allele"]

    # Building our custom hash:
    variation_data = {
        "reference": ref,
        "rsID" : variationData["name"],
        "minorAllele" : variationData["minor_allele"],
        "ancestralAllele" : variationData["ancestral_allele"],
        "synonyms" : variationData["synonyms"],
        "MAF" : variationData["MAF"],
        "varClass" : variationData["var_class"],
        "consequence" : variationData["most_severe_consequence"]
    }

    # Extracting frequencies:
    pop_freq = {}
    for pops in variationData["populations"]:
        pop_name = pops["population"].split(":")
        try:
            if pop_name[0] == "1000GENOMES" and pop_name[1] == "phase_3":
                try:
                    pop_freq[pop_name[2]][pops["allele"]] = pops["frequency"]
                except:
                    pop_freq[pop_name[2]] = {pops["allele"] : pops["frequency"]}
        except:
            continue
    # If not 1000 genome frequencies are available for the given variation:
    if len(pop_freq) == 0:
        pop_freq = "Variation was not found in the 1000 Genomes data."

    # Extracting mapping:
    for mapping in variationData["mappings"]:
        if len(mapping["seq_region_name"]) > 2:
            continue
        else:
            variation_data["chromosome"] = mapping["seq_region_name"]
            variation_data["start"] = mapping["start"]
            variation_data["end"] = mapping["end"]
            variation_data["allele_string"] = mapping["allele_string"]
            variation_data["location"] = mapping["location"]

    # 'alternativeAllele': alt,
    if variationData["minor_allele"] == ref:
        all_alleles = variation_data["allele_string"].split("/")
        variation_data["alternativeAllele"] = all_alleles[1]
    else:
        variation_data["alternativeAllele"] = variationData["minor_allele"]

    # Extracting phenotype information:
    if len(variationData["phenotypes"]) > 0: variation_data["phenotypes"] = []
    for pt in variationData["phenotypes"]:
        trait = pt["trait"] if "trait" in pt else "-"
        source = pt["source"] if "source" in pt else "-"
        allele = pt["risk_allele"] if "risk_allele" in pt else "-"

        if "risk_allele" in pt:
            trait += "("+pt["risk_allele"]+")"

        variation_data["phenotypes"].append({
            "trait": trait,
            "source": source,
            })



    # Generating SNP ID:
    variation_data["SNP_ID"] = "chr" + str(variation_data["chromosome"]) + ":" + str(variation_data["start"])

    # Extracting clinical significance:
    try:
        for cs in variationData["clinical_significance"]:
            if cs != "other" and cs != "not provided":
                try:
                    variation_data["clinical_significance"].append(cs)
                except:
                    variation_data["clinical_significance"] = [cs]
    except:
        pass


    # Now we have to generate the gerp and gwava scores.
    # variation_data = get_GWAVA_score(variation_data)

    return (pop_freq, variation_data)

# Function to retrieve Exome Aggregation Consortium (ExAC) allele frequencies:
def get_ExAC_frequencies(variant_data):
    '''
    This function runs a tabix query to retrieve allele frequency data published in the Exome Aggreagtion Consortium.
    The query is based on the variation start and end coordinates, but also involves an allele checking to make sure
    We don't retrieve overlapping but non identical variations.

    The returned dictionary contains allele frequencies and counts for all reported populations for all alternative alleles.

    It is important that the alternative alleles are reported even if the minor allele is reference allele! For this
    reason we test the alleles in the get_variant_info() function.

    '''
    # Required parameters:
    queried_var = variant_data["location"]
    start = variant_data["start"]
    end = variant_data["end"]
    reference  = variant_data['reference']

    # In the Ensembl annotation, sometimes, the alternative allele is set as NoneType.
    # In such cases the alternative has to be extracted:
    alternative = variant_data['alternativeAllele']
    if alternative is None:
        for allele in variant_data['allele_string'].split("/"):
            if allele != reference:
                alternative = allele
                break

    # ExAC file is stored in the config.py file:
    ExAC_file = config.EXAC_FILE

    # Check if exac file exists or not:
    if not os.path.isfile(ExAC_file):
        return "ExAC datafiles were not found."

    query = "bash -O extglob -c \'/software/hgi/pkglocal/tabix-git-1ae158a/bin/tabix %s %s\'" %(ExAC_file, queried_var)

    # Submit query:
    output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Initializing returned variable:
    Exac_parsed = {}
    for line in output.stdout.readlines():
        fields = line.strip().split("\t")
        # Checking if the line is indeed the requested variant:
        if (str(fields[1]) == str(start) or str(fields[1]) == str(end)) and (fields[3] == reference)  and (alternative in fields[4]):
            ExAC_data = {}
            alternative_alleles = fields[4].split(",")

            # Reading data from line:
            for detail in fields[7].split(";"):
                try:
                    (key, value) = detail.split("=")
                    ExAC_data[key] = value.split(",")
                except:
                    pass

            # Parsing data:
            ExAC_pops = ["AFR", "ALL", "FIN", "AMR", "EAS", "SAS", "OTH", "NFE"]
            for i in range(len(alternative_alleles)):
                Exac_parsed[alternative_alleles[i]] = {}

                # Now, we have to extract all allele counts and frequencies:
                for population in ExAC_pops:
                    if population == "ALL": pop = ""
                    else: pop = "_" + population

                    freqency = float(ExAC_data["AC"+pop][i]) / int(ExAC_data["AN"+pop][0])
                    Exac_parsed[alternative_alleles[i]][population] = {"count": ExAC_data["AC"+pop][i], "frequency":freqency}
        else:
            pass
    if len(Exac_parsed) > 0:
        return Exac_parsed
    else:
        return "Variation could not be found in the ExAC dataset."

# This function checks if rsID is exists for a given snp id
def get_rsID(SNP_ID):
    '''
    This is the first step of the variation annotation, when we test if the provided SNP ID is
    proper, and if it overlaps with an already known rsID.

    If everything is OK, we return with the population and the variation data
    '''

    # If rsID is provided than we don't do all the crap:
    if "rs" in SNP_ID:
        return get_variant_info(SNP_ID, "NA")

    # Step1: parse input data: assign chromosome and position and the allele string:
    coordinate = SNP_ID.split("_")[0]
    alleles = SNP_ID.split("_")[1]

    chromosome = coordinate.split(":")[0]
    if "chr" in chromosome: chromosome = chromosome[3:]
    position = coordinate.split(":")[1]

    allele1 = alleles.split("/")[0]
    allele2 = alleles.split("/")[1]

    # Step2: checking reference sequence:
    URL = config.REST_URL + "/sequence/region/human/%s:%s..%s:1?content-type=text/plain" % (chromosome, position, position)

    base = submit_REST(URL)["seq"]

    if allele1 == base:
        ref = allele1
        alt = allele2
    elif allele2 == base:
        ref = allele2
        alt = allele1
    else:
        print "[Error] There is a problem with the input: %s, none of the alleles are matching with the reference allele (%s).\n" %(SNP_ID, base)

    # Checking for overlapping variations:
    URL = config.REST_URL + "/overlap/region/human/%s:%s-%s?feature=variation" % (chromosome, position, position)
    variations = submit_REST(URL)

    rsID = ""
    for variation in variations:
        try:
            if ref in variation["alt_alleles"] and alt in variation["alt_alleles"]:
                rsID = variation["id"]
        except:
            if ref in variation["alleles"] and alt in variation["alleles"]:
                rsID = variation["id"]

    variant_data = []
    if "rs" in rsID:
        return get_variant_info(rsID, ref)
    else:
        # checking variation class:
        if len(ref) == len(alt):
            var_class = "SNP"
        elif len(ref) > len(alt):
            var_class = "DEL"
        elif len(ref) < len(alt):
            var_class = "INS"

        return [{},{
             'reference' : ref,
             'MAF': '-',
             'SNP_ID': coordinate,
             'allele_string': alleles,
             'ancestralAllele': ref,
             'chromosome': chromosome,
             'consequence': '',
             'end': int(position),
             'location': coordinate+'-'+position,
             'minorAllele': alt,
             'alternativeAllele': alt,
             'rsID': "-",
             'start': int(position),
             'synonyms': [],
             'varClass': var_class

        }]

# A function to retrieve publications where the given variation was mentioned:
def get_pubmedIDs(rsID):
    '''
    This function returns the list of pubmed ID-s of those publications where the given rsID was mentioned.

    Up to 1000 IDs are returned. If a variation was mentioned in more papers, some won't show up in the table.
    '''

    # Formulation of the query syntax:
    URL = config.PUBMED_URL_VAR % (rsID)

    # Submitting request, and parsing result:
    r = requests.get(URL)
    decoded = r.json()

    # Returning the list of IDs:
    pubmed_IDs = decoded["esearchresult"]["idlist"]

    # Before we return the list, let's return all the titles for all articles:
    publication_data = {}
    for ID in pubmed_IDs:

        URL = config.PUBMED_URL_PMID % (ID)
        r = requests.get(URL)
        decoded = r.json()

        # Extracting data:
        publication_data[ID] = {
            "firstAuthor" : decoded["result"][ID]['sortfirstauthor'],
            "title" : decoded["result"][ID]['title'],
            "journal" : decoded["result"][ID]['fulljournalname'],
            "year" : decoded["result"][ID]['epubdate'].split(" ")[0],
            "URL" : "http://www.ncbi.nlm.nih.gov/pubmed/"+ID,
        }


    return publication_data

# Checking overlapping regulatory elements:
def get_Ensembl_regulation(variant_data):
    '''
    Based on chromosome and position, this function downloads
    '''
    # We slightly extend our search for +/- 10bp:
    start = variant_data["start"] - 10
    end = variant_data["end"] + 10
    chromosome = variant_data["chromosome"]

    URL = config.REST_URL + "/overlap/region/human/%s:%s-%s?feature=regulatory;content-type=application/json" % (chromosome, start, end)
    content = submit_REST(URL)

    regulatory_table = {}

    # Reading data from the returned variable:
    for regulatory in content:

        # Calculating distance:
        distance = min(
            abs(start - regulatory["start"]),
            abs(start - regulatory["end"]),
            abs(end - regulatory["start"]),
            abs(end - regulatory["end"])
        )
        if variant_data["start"] > regulatory["start"] and variant_data["start"] < regulatory["end"]: distance = 0

        # Adding data to dict:
        regulatory_table[regulatory["ID"]] = {
            "ID" : regulatory["ID"],
            "cell_type": regulatory["cell_type"],
            "description": regulatory["description"],
            "distance": distance
        }

    return regulatory_table

# Based on a genomic location, this function retrieves a list of genes within a 1Mbp window.
def get_gene_list(variant_data):
    '''
    Based on the submitted chromosome and position, this function returns
    all overlapping genes within 1Mbp window using the REST API of Ensembl

    The returned dictionary contains the basic information
    '''

    # Creating a 1Mbp window around the variation:
    chromosome = variant_data["chromosome"]
    position = variant_data["start"]

    end = position + config.WINDOW
    start = position - config.WINDOW
    if start < 0: start = 0

    # Now, find overlapping genes:
    URL = config.REST_URL + "/overlap/region/human/%s:%s-%s?feature=gene;content-type=application/json" %(chromosome, start, end)

    overlapping_genes = submit_REST(URL)

    gene_list = {}
    for gene in overlapping_genes:

        # This dictionary will be just added into the returned bid dataset:
        gene_details = {
            "start" : gene["start"],
            "end" : gene["end"],
            "strand" : gene["strand"],
            "name" : gene["external_name"],
            "description" : gene["description"],
            "biotype" : gene["biotype"],
            "ID" : gene["id"]
        }

        # Calculate distance:
        d_from_start = abs(position - gene["start"])
        d_from_end = abs(position - gene["end"])

        distance = min(d_from_start, d_from_end)
        if position >= gene["start"] and position <= gene["end"]:
            distance = 0

        if distance not in gene_list.keys():
            gene_list[distance] = []

        gene_details["distance"] = distance

        # calculate orientation:
        gene_details["orientation"] = "upsteram"
        if gene_details["distance"] == 0:
            gene_details["orientation"] = "overlapping"
        elif gene["strand"] == 1 and d_from_end < d_from_start:
            gene_details["orientation"] = "downstream"
        elif gene["strand"] == -1 and d_from_end > d_from_start:
            gene_details["orientation"] = "downstream"

        # Adding novel data into the dict:
        gene_list[distance].append(gene_details)

    return gene_list

# This function is for running a regulatory rings:
def get_regulation(variant_data):
    '''
    This function returns all overlapping regulatory features within 2kb of the variation.
    The function is based on a local bedtool query on a datafile downloaded from Ensembl
    (data mainly sourced from Encode/Epigenomics Roadmap data)

    Input: chromosome/position data is read from the variant_data variable previously created.

    Putput: dictionary with all regulatory feature classes as keys, and types, and cell types as
    values.

    '''

    # Cell type abbreviations are retrieved from Ensembl, which are replaced
    # by conventional names for better readibility.
    # Not all cell types are included from tier-3 cell types. <- should be updated
    # if needed. Cell types are moved to config.py
    CellTypeFinder = config.CellTypeFinder

    chromosome = variant_data["chromosome"]
    start = variant_data["start"] - 500
    end = variant_data["end"] + 500

    # folders of regulatory files:
    regulatoryFile = config.REGULATORY_FILE

    # Seraching for modification:
    query = "bash -O extglob -c \'/nfs/team144/software/bedtools2/bin/intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\n\") -b %s -sorted\'" %  (
        chromosome, start, end, regulatoryFile)

    # Submit query:
    output = subprocess.Popen(query.strip(),
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)

    # Once the queries are returned, we have to parse the output:
    annotations = []
    for line in output.stdout.readlines():
        x = {}
        for field in line.strip().split("\t")[6].split(";"):
            x[field.split("=")[0]] = field.split("=")[1]
        annotations.append(x)

    # Parsing returned information:
    parsed_annotations = {}
    for annotation in annotations:

        # Parsing modification type:
        if not annotation["Class"] in parsed_annotations.keys():
            parsed_annotations[annotation["Class"]] = {}
        if not annotation["Name"] in parsed_annotations[annotation["Class"]].keys():
            parsed_annotations[annotation["Class"]][annotation["Name"]] = []

        # Parsing cell type:
        try:
            cell = annotation["Cell_type"]
        except:
            cell = annotation["Alias"].split("_")[0]
        try:
            cell = CellTypeFinder[cell]
        except:
            pass

        # We just adding the cell type:
        if not cell in parsed_annotations[annotation["Class"]][annotation["Name"]]:
            parsed_annotations[annotation["Class"]][annotation["Name"]].append(cell)

    return (parsed_annotations)

# Retrieve UK10K allele frequency:
def get_UK10K_frequencies(variant_data):
    '''
    This function queries the uk10k public data for allele frequency and allele counts
    It need the previously generate variant_data dictionary with all details:
        - chromosome, variant location, start, end alternative, reference, allele string, rsID

    Returned data: dictionary with the following structure:
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
    In the genotype counts, the numbers in parenthese indicate the expected
    counts assuming perfect HW eq

    If error occurred: string explaining the situation
        1) Missing datafiles.
        2) Matching variations were not found.

    '''

    chromosome = variant_data["chromosome"]
    queried_var = variant_data["location"]
    start = variant_data["start"]
    end = variant_data["end"]
    alternative = variant_data['alternativeAllele']
    reference  = variant_data["reference"]
    allele_string = variant_data["allele_string"]
    rsID = variant_data["end"]

    # Check if datafile is exists, and return 0 if not:
    return "UK10K datafiles were not found."

    # tabix indexed vcf file, in which we are looking for the variation:
    UK10K_vcf = config.UK10K_FILE % (chromosome)

    # Check if datafile is exists, and return 0 if not:
    if not os.path.isfile(UK10K_vcf):
        return "UK10K datafiles were not found."

    # submit bcftools query:
    query = "bash -O extglob -c \'/software/hgi/pkglocal/bcftools-1.2/bin/bcftools query -f \"%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n\" -r %s %s\'" %(queried_var, UK10K_vcf)
    output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Initializing returned variable:
    genotype_count = {}
    for line in output.stdout.readlines():

        # from each lines, the info fields are extracted for checking:
        fields = line.strip().split("\t")
        out_chromosome = fields.pop(0)
        out_position = fields.pop(0)
        out_ID = fields.pop(0)
        out_reference = fields.pop(0)
        out_alternative = fields.pop(0)

        # Now we test if the variation was found in the UK10K database:
        if rsID == out_ID or (chromosome == out_chromosome and
            (int(out_position) == start or int(out_position) == end) and
            out_alternative in allele_string and out_reference in allele_string):

            # Now we parse the genotypes:
            genotype_count = {
                "Allele_counts": {
                    out_reference : 0,
                    out_alternative : 0
                },
                "Genotype_counts" : {
                    out_reference+out_reference : 0,
                    out_alternative+out_reference : 0,
                    out_alternative+out_alternative : 0,
                },
                "Allele_frequencies": {
                    out_reference : 0,
                    out_alternative : 0
                }
            }

            n = 0 # Counting samples in row
            for sample in fields:
                try:
                    gt = sample.split("=")[1]
                    (a1, a2) = gt.split("|")

                    n += 1

                    # Parsing genotypes:
                    if a1 == "1" and a2 == "1":
                        genotype_count["Allele_counts"][out_alternative] += 2
                        genotype_count["Genotype_counts"][out_alternative+out_alternative] += 1
                    if a1 == "0" and a2 == "0":
                        genotype_count["Allele_counts"][out_reference] += 2
                        genotype_count["Genotype_counts"][out_reference+out_reference] += 1
                    else:
                        genotype_count["Allele_counts"][out_reference] += 1
                        genotype_count["Allele_counts"][out_alternative] += 1
                        genotype_count["Genotype_counts"][out_alternative+out_reference] += 1

                except:
                    pass # Unexpected field. Skip

            # Calculating the allele frequencies:
            genotype_count["Allele_frequencies"][out_reference] = round(float(genotype_count["Allele_counts"][out_reference]) /(
                                    genotype_count["Allele_counts"][out_alternative] + genotype_count["Allele_counts"][out_reference]), 3)
            genotype_count["Allele_frequencies"][out_alternative] = round(float(genotype_count["Allele_counts"][out_alternative]) / (
                                    genotype_count["Allele_counts"][out_alternative] + genotype_count["Allele_counts"][out_reference]), 3)

            # Calculating expected number of genotypes (based on pop size and HW):
            exp_AA = (genotype_count["Allele_frequencies"][out_reference] ** 2) * n
            exp_aa = (genotype_count["Allele_frequencies"][out_alternative] ** 2) * n
            exp_Aa = 2 * (genotype_count["Allele_frequencies"][out_alternative] * genotype_count["Allele_frequencies"][out_reference]) * n

            # Updating genotype counts with the expected values:
            genotype_count["Genotype_counts"][out_alternative+out_alternative] = '%s (%s)'  %(genotype_count["Genotype_counts"][out_alternative+out_alternative],
                                                                                           round(exp_aa, 2))
            genotype_count["Genotype_counts"][out_reference+out_reference] = '%s (%s)'  %(genotype_count["Genotype_counts"][out_reference+out_reference],
                                                                                          round(exp_AA, 2))
            genotype_count["Genotype_counts"][out_alternative+out_reference]= '%s (%s)'  %( genotype_count["Genotype_counts"][out_alternative+out_reference],
                                                                                          round(exp_Aa, 2))

            # If we found only one matching line, we return data:
            return genotype_count

        else:
            pass # Going to the next line.


    # Once we reached the end of the lines without finding the correct line, then we
    return "Variation was not found in the UK10K databaset."

# Found eQTL for a given variation:
def get_GTEx_genes(variant_data):
    '''
    This function was written to check if a variation is an eQTL of any genes.
    Input: genomic coordinates stored in the previously generated variant data dictionary
    Output: dictionary with tissues as keys, then every gene in that tissue is an element of an array
    If a problem occurred, a the appropriate error message returned as string.
    '''

    # GTEx file is listed in the config.py file:
    GTEx_file = config.GTEX_FILE_VAR

    # Check if datafile is exists, and return 0 if not:
    if not os.path.isfile(GTEx_file):
        return "GTEx datafile was not found in this location: %s" % GTEx_file

    # Extract data for bedtools query:
    chromosome = variant_data["chromosome"]
    start = variant_data["start"]
    end = variant_data["end"]

    # submit bcftools query:
    query = "bash -O extglob -c \'/software/hgi/pkglocal/tabix-git-1ae158a/bin/tabix %s chr%s:%s-%s\'" %(GTEx_file, chromosome, start, end + 1)
    output = subprocess.Popen(query.strip(), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    GTEx_data = {}
    for line in output.stdout.readlines():
        fields = line.strip().split("\t")
        if fields[3] not in GTEx_data:
            GTEx_data[fields[3]] = []

        if fields[4] == variant_data["rsID"] or fields[4] in variant_data["synonyms"]:
            try:
                parts = str(float(fields[11])).split("e")
                pval = parts[0][0:5]+"e"+parts[1]
            except:
                pval = round(float(fields[11]), 3)
            beta = round(float(fields[9]), 3)
            sd = round(float(fields[10]), 3)
            tissue = fields[3].replace("_", " ")
            GTEx_data[fields[3]].append({
                "rsID" : fields[4],
                "tissue" : tissue,
                "gene_name" : fields[5],
                "gene_ID" : fields[6],
                "distance" : fields[8],
                "beta" : beta,
                "sd" : sd,
                "biotype" : fields[7],
                "pvalue"  : pval
            })

    # Removing empty keys of the dictionary:
    for tissue in GTEx_data.keys():
        if len(GTEx_data[tissue]) == 0: GTEx_data.pop(tissue)

    if len (GTEx_data) == 0:
        GTEx_data = "The queried variation was not found in the GTEx dataset."

    return GTEx_data
