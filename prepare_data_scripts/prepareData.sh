#!/bin/bash

function usage {
    echo ""
    echo "Usage:" $(basename $0) "-o <output dir> { -h -g -x -r}"
    echo " -h : this help message"
    echo " -g -x -r are optional, for selecting which data to prepare (GWAS, GTEx or Regulation)"
}

OPTIND=1

prepgwas=0
prepgtex=0
prepreg=0

while getopts "gxro:h" optname; do
    case "$optname" in
        "g" ) prepgwas=1;;
        "x" ) prepgtex=1;;
        "r" ) prepreg=1;;
        "o" ) out="${OPTARG}";;
        "h" ) usage ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

if [ "$prepgwas" -eq 0 ] && [ "$prepgtex" -eq 0 ] && [ "$prepreg" -eq 0 ];then
    prepgwas=1
    prepgtex=1
    prepreg=1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CDIR=$(pwd)

reg_script="$DIR"/"prepareRegulation.sh"
gwas_script="$DIR"/"prepareGWAS.py"
gtex_script="$DIR"/"prepareGTEx.py"

regdir=${regdir/%\/}
out=${out/%\/}

# ------------------------------- URLs ------------------------------------

gtexURL="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar"
gtexlookupURL="https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
gtexgencodeURL="https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"
gwasURL="https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
ensregURL="ftp://ftp.ensembl.org/pub/release-100/regulation/homo_sapiens/RegulatoryFeatureActivity/"

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating output directories"
mkdir -p "$out/gwas"
mkdir -p "$out/gtex"
mkdir -p "$out/regulation"
mkdir -p "$out/temp"

# ----------------------- DOWNLOADING DATA --------------------------------

if [[ "$prepgtex" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx data"
    wget --quiet -c "$gtexURL" -O "$out/temp/GTEx.tar"

    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx lookup data"
    wget --quiet -c "$gtexlookupURL" -O "$out/gtex/GTEx.lookup.txt.gz"

    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx Gencode data"
    wget --quiet -c "$gtexgencodeURL" -O "$out/gtex/gencode.gtf"
    gzip -f "$out/gtex/gencode.gtf"
fi

if [[ "$prepgwas" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GWAS catalog"
    wget --quiet -c "$gwasURL" -O "$out/gwas/gwas_full.tsv"
    gzip -f "$out/gwas/gwas_full.tsv"
fi

if [[ "$prepreg" -eq 1 ]];then
    cd "$out/temp"
    mkdir -p regulation
    cd regulation
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading Ensembl Regulation"
    lftp -c "open $ensregURL; mirror -P . ." > /dev/null 2> /dev/null
    cd "$CDIR"
fi

# -------------------------- REGULATION -----------------------------------

if [[ "$prepreg" -eq 1 ]];then
    regdir="$out/temp/regulation"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating regulation file"
    "$reg_script" -i "$regdir" -o "$out/regulation/regulation.bed"
fi

# ----------------------------- GTEx --------------------------------------

if [[ "$prepgtex" -eq 1 ]];then
    gtex="$out/temp/GTEx.tar"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GTEx file"
    "$gtex_script" -i "$gtex" -g "$out/gtex/gencode.gtf.gz" -l "$out/gtex/GTEx.lookup.txt.gz" | sort -k1,1 -k2,2n > "$out/gtex/gtex.bed"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Compressing GTEx file"
    bgzip -f "$out/gtex/gtex.bed"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Indexing GTEx file"
    tabix -f -p bed "$out/gtex/gtex.bed.gz"
fi

# ----------------------------- GWAS ---------------------------------------

if [[ "$prepgwas" -eq 1 ]];then
    gwas="$out/gwas/gwas_full.tsv.gz"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GWAS file"
    PYTHONPATH=$(dirname "$DIR") "$gwas_script" -i "$gwas" | gzip - > "$out/gwas/gwas.tsv.gz"
fi

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Deleting temporary data"
rm -rf "$out/temp/*"

exit 0
