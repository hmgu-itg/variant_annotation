# batch sizes for POST queries
VARIANT_RECODER_POST_MAX=200
VARIATION_POST_MAX=200
VEP_POST_MAX=200

# window around a variant overlapping ENSEMBL regulatory elements (regulation::getRegulation())
REG_WINDOW=2000
# window around a variant to get variants annotated with phenotypes from (variant::getVariantsWithPhenotypes())
PHENO_WINDOW=500000
# window around a genomic position to get a list of genes overlapping it (gene::getGeneList())
GENE_WINDOW=1000000
# window around a variant to get GWAS signals from
GWAS_WINDOW=500000

GNOMAD_URL="https://gnomad.broadinstitute.org/api"
#UNIPROT_URL="http://www.uniprot.org/uniprot/"
UNIPROT_URL="https://rest.uniprot.org/uniprotkb/"

PUBMED_URL_VAR="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%%22%s%%22%%5BAll%%20Fields%%5D&retmode=json&retmax=1000"
PUBMED_URL_PMID="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&rettype=abstract&id=%s"
PUBMED_URL="https://pubmed.ncbi.nlm.nih.gov/"

# links for various phenotype sources 
CLINVAR_URL="https://www.ncbi.nlm.nih.gov/clinvar/?term="
NHGRI_URL="https://www.ebi.ac.uk/gwas/search?query="
ENSEMBL_PHENO_URL="http://ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;v=%s;vdb=variation"

# prepared local files
REGULATORY_FILE_SFX="regulation/regulation.bed.gz"
GWAS_FILE_VAR_SFX="gwas/gwas.tsv.gz"
GWAS_FILE_SFX="gwas/gwas_full.tsv.gz"
GTEX_BED_SFX="gtex/gtex.bed.gz"
GXA_TSV_SFX="gxa/GXA.tsv.gz"

# VEP consequences severity
VEP_CONSEQUENCES={
    "transcript_ablation": 39,
    "splice_acceptor_variant": 38,
    "splice_donor_variant": 37,
    "stop_gained": 36,
    "frameshift_variant": 35,
    "stop_lost": 34,
    "start_lost": 33,
    "transcript_amplification": 32,
    "inframe_insertion": 31,
    "inframe_deletion": 30,
    "missense_variant": 29,
    "protein_altering_variant": 28,
    "splice_region_variant": 27,
    "splice_donor_5th_base_variant": 26,
    "splice_donor_region_variant": 25,
    "splice_polypyrimidine_tract_variant": 24,
    "incomplete_terminal_codon_variant": 23,
    "start_retained_variant": 22,
    "stop_retained_variant": 21,
    "synonymous_variant": 20,
    "coding_sequence_variant": 19,
    "mature_miRNA_variant": 18,
    "5_prime_UTR_variant": 17,
    "3_prime_UTR_variant": 16,
    "non_coding_transcript_exon_variant": 15,
    "intron_variant": 14,
    "NMD_transcript_variant": 13,
    "non_coding_transcript_variant": 12,
    "upstream_gene_variant": 11,
    "downstream_gene_variant": 10,
    "TFBS_ablation": 9,
    "TFBS_amplification": 8,
    "TF_binding_site_variant": 7,
    "regulatory_region_ablation": 6,
    "regulatory_region_amplification": 5,
    "feature_elongation": 4,
    "regulatory_region_variant": 3,
    "feature_truncation": 2,
    "intergenic_variant": 1
    # "transcript_ablation"      : 36,
    # "splice_acceptor_variant"  : 35,
    # "splice_donor_variant"     : 34,
    # "stop_gained"              : 33,
    # "frameshift_variant"       : 32,
    # "stop_lost"                : 31,
    # "start_lost"               : 30,
    # "transcript_amplification" : 29,
    # "inframe_insertion"        : 28,
    # "inframe_deletion"         : 27,
    # "missense_variant": 26,
    # "protein_altering_variant": 25,
    # "splice_region_variant": 24,
    # "incomplete_terminal_codon_variant": 23,
    # "start_retained_variant": 22,
    # "stop_retained_variant": 21,
    # "synonymous_variant": 20,
    # "coding_sequence_variant": 19,
    # "mature_miRNA_variant": 18,
    # "5_prime_UTR_variant": 17,
    # "3_prime_UTR_variant": 16,
    # "non_coding_transcript_exon_variant": 15,
    # "intron_variant": 14,
    # "NMD_transcript_variant": 13,
    # "non_coding_transcript_variant": 12,
    # "upstream_gene_variant": 11,
    # "downstream_gene_variant": 10,
    # "TFBS_ablation": 9,
    # "TFBS_amplification": 8,
    # "TF_binding_site_variant": 7,
    # "regulatory_region_ablation": 6,
    # "regulatory_region_amplification": 5,
    # "feature_elongation": 4,
    # "regulatory_region_variant": 3,
    # "feature_truncation": 2,
    # "intergenic_variant": 1
}

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
