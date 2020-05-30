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

##
## Shared parameters
##
REST_URL    = 'http://grch37.rest.ensembl.org'
APPRIS_URL  = 'http://apprisws.bioinfo.cnio.es:80/rest/exporter/id/homo_sapiens/%s?format=json&db=hg19'
WINDOW      = 500000

##
## Templates:
##
GENE_TEMPLATE   = '/nfs/team144/ds26/FunctionalAnnotation/v2.2/template_gene.html'
VAR_TEMPLATE    = '/nfs/team144/ds26/FunctionalAnnotation/v2.2/template.html'

##
## Parameters for genes:
##
OMIM_URL        = 'http://api.omim.org/api/entry?'
OMIM_APIKEY     = 'apiKey=MWkzg-DYShWyrfrbqUHL4g&format=python&' # will expire on March 21st, 2017.
GWAS_FILE       = '/nfs/team144/ds26/FunctionalAnnotation/GWAS_positive_controls/full'
UNIPROT_URL     = 'http://www.uniprot.org/uniprot/'
GXA_URL         = 'http://www.ebi.ac.uk/gxa/widgets/heatmap/multiExperiment.tsv?geneQuery=%s&propertyType=bioentity_identifier&species=homo%%20sapiens'
GTEX_FILE_GENES = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/GTEx/GTEx_for_genes.bed.gz'

##
## Parameters for variations:
##
GWAVA_DIR       = "/lustre/scratch113/teams/zeggini/users/ds26/GWAVA/gwava_release"
GWAS_FILE_VAR   = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/GWAS_catalog/gwas_catalog_2016.02.10.bed.gz'
EXAC_FILE       = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/ExAC_rel0.3/ExAC.r0.3.nonpsych.sites.vcf.gz'
PUBMED_URL_VAR  = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%%22%s%%22%%5BAll%%20Fields%%5D&retmode=json&retmax=1000'
PUBMED_URL_PMID = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract&id=%s'
REGULATORY_FILE = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/regulation/Ensembl_regulation_2015.11.13_AnnotatedFeatures.bed.gz'
REG_FEAT_FILE   = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/regulation/ALL_Ensembl_regulatory_features_resorted.bed.gz'
UK10K_FILE      = '/nfs/vertres18/projects/uk10k/RELEASE/UK10K_COHORT/REL-2012-06-02/v3/%s.beagle.anno.csq.shapeit.20131101.vcf.gz'
UK10K_VERSION   = 'UK10K data release 2013 Nov. 11, version: v3.'
GTEX_FILE_VAR   = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/GTEx/GTEx_all_tissue.bed.gz'

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

##
## Variables used for parsing outputs:
##
CellTypeFinder = {
    "K562" : "Myelogenous leukemia line",
    "H1-hESC" : "Embryonic stem cells",
    "GM12878" : "Lymphoblastoid cell line",
    "A549" : "Adenocarcinomic alveolar epithelial cells",
    "CD20+" : "B-cells CD20+",
    "CD20+_RO01778" : "B-cells CD20+ (RO01778)",
    "CD20+_RO01794" : "B-cells CD20+ (RO01794)",
    "H1-neurons" : "Neuronal cell line",
    "HeLa-S3" : "Cervical cancer",
    "HepG2" : "Liver hepatocellular carcinoma",
    "HUVEC" : "Umbilical vein/vascular endothelium",
    "IMR90" : "Fibroblast",
    "LHCN-M2" : "Skeletal myoblasts",
    "MCF-7" : "Breast cancer cell line",
    "Monocytes-CD14+" : "Monocytes-CD14+",
    "Monocytes-CD14+_RO01746" : "Monocytes-CD14+ (RO01746)",
    "Monocytes-CD14+_RO01826" : "Monocytes-CD14+ (RO01826)",
    "SK-N-SH" : "Neuroblastoma cell line",
    "HSMMtube" : "Skeletal muscle myotubes",
    "H1ESC" : "Embryonic stem cells",
    "DND-41" : "T cell leukemia",
    "HMEC" : "Mammary epithelial cells",
    "Osteobl" : "Osteoblast",
    "HSMM" : "Skeletal muscle myoblasts",
    "NH-A" : "Astrocytes",
    "NHDF-AD" : "Dermal fibroblasts",
    "NHEK" : "Epidermal keratinocytes",
    "NHLF" : "Lung fibroblasts",
}

population_names = {
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
    "AMR" : "Total Americas Ancestry",
    "ALL" : "Total 1000 Genomes dataset"
}
