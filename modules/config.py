'''
Parameter file for the Variant annotation tool. Parameters are strored in a dictionary
that can be accessed from any modules that import this one.

Version: 2.2.0 Last modified: 2016.02.24

Usage:

    * This file is easy to extend, just add new key/values separated by
      equal sign. Values can be access from any function by referring to its
      key: value = GlobalParameters['key']
    * Comment lines starting with # character will be ignored. Empty lines
      will also be ignored.
    * Time to time the OMIM API key have to be re-newed.
    * Keys should not be, values should be double or single quoted.

'''

##
## Version
##
VERSION = "2.2.0"

## batch size for POST queries
BATCHSIZE=100

# window around a variant overlapping ENSEMBL regulatory elements (regulation::getRegulation())
REG_WINDOW=2000
# window around a variant to get variants annotated with phenotypes from (variant::getVariantsWithPhenotypes())
PHENO_WINDOW=500000
# window around a genomic position to get a list of genes overlapping it (gene::getGeneList())
GENE_WINDOW=1000000
# window around a variant to get GWAS signals from
GWAS_WINDOW=500000

GNOMAD_URL="https://gnomad.broadinstitute.org/api"
APPRIS_URL="http://apprisws.bioinfo.cnio.es:80/rest/exporter/id/homo_sapiens/%s?format=json&db=hg19"

OMIM_URL        = 'http://api.omim.org/api/entry?'
OMIM_APIKEY     = 'apiKey=MWkzg-DYShWyrfrbqUHL4g&format=python&' # expired
UNIPROT_URL     = 'http://www.uniprot.org/uniprot/'
GXA_URL_PREFIX  = "https://www.ebi.ac.uk/gxa/genes/"
GXA_URL_SUFFIX  = "?bs=%7B%22homo%20sapiens%22%3A%5B%22ORGANISM_PART%22%5D%7D&ds=%7B%22kingdom%22%3A%5B%22animals%22%5D%7D#baseline"
PUBMED_URL_VAR  = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%%22%s%%22%%5BAll%%20Fields%%5D&retmode=json&retmax=1000'
PUBMED_URL_PMID = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract&id=%s'

##
## Anchor tags for Ensembl and other sources:
##
ENSEMBL_VAR     = '<a href="http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=%s;vdb=variation">%s</a>'
ENSEMBL_GENE    = '<a href="http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=%s">%s</a>'
ENSEMBL_REGION  = '<a href="http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;r=%s">%s</a>'
CLINVAR_LINK    = '<a href="http://www.ncbi.nlm.nih.gov/clinvar/?database=clinvar&term=%s">%s</a>'
OMIM_LINK       = '<a href="http://www.omim.org/search?index=entry&start=1&limit=10&search=%s">%s</a>'
GWAS_CAT_LINK   = '<a href="https://www.ebi.ac.uk/gwas/search?query=%s">%s</a>'
PUBMED_LINK     = '<a href="http://www.ncbi.nlm.nih.gov/pubmed/%s">%s</a>'
GTEX_LINK       = '<a href="http://www.gtexportal.org/home/snp/%s">GTEx site</a>'

GnomadPopulations={
    "AFR": "African American",
    "AMI": "Amish",
    "AMR": "Latino",
    "ASJ": "Ashkenazi Jews",
    "EAS": "East Asian",
    "FIN": "Finnish",
    "NFE": "Non-Finnish European",
    "OTH": "Other",
    "SAS": "South Asian"
}

PopulationNames = {
    "CDX" : "Chinese Dai in Xishuangbanna",
    "CHB" : "Han Chinese in Bejing",
    "JPT" : "Japanese in Tokyo, Japan",
    "KHV" : "Kinh in Ho Chi Minh City",
    "CHS" : "Southern Han Chinese",
    "EAS" : "Total East Asian Ancestry",
    "BEB" : "Bengali",
    "GIH" : "Gujarati Indian",
    "ITU" : "Indian Telugu",
    "PJL" : "Punjabi in Lahore",
    "STU" : "Sri Lankan Tamil",
    "SAS" : "Total South Asian Ancestry",
    "ASW" : "African American",
    "ACB" : "African Caribbean",
    "ESN" : "Esan",
    "MAG" : "Gambian",
    "LWK" : "Luhya",
    "MSL" : "Mende",
    "YRI" : "Yoruba",
    "AFR" : "Total African Ancestry",
    "GBR" : "British",
    "FIN" : "Finnish",
    "IBS" : "Iberian",
    "TSI" : "Toscani",
    "CEU" : "Western European ancestry",
    "EUR" : "Total European Ancestry",
    "CLM" : "Colombian",
    "MXL" : "Mexican Ancestry",
    "PEL" : "Peruvian",
    "PUR" : "Puerto Rican",
    "AMR" : "Total American Ancestry",
    "ALL" : "Total 1000 Genomes dataset"
}

# how many tissues with highest expression values to report from GXA
GXA_HIGHEST=10
