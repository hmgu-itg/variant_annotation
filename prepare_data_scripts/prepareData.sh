#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -o <output dir>"
#    echo "Usage: $0 -x <input GTEx.tar>"
#    echo "          -o <output dir>"
}

OPTIND=1
while getopts "o:" optname; do
    case "$optname" in
#        "x" ) gtex="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

# if [[ ! -f "$gtex" ]];then
#     echo $(date '+%d/%m/%Y %H:%M:%S') "ERROR: GTEx file ($gtex) does not exist"
#     exit 1
# fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CDIR=$(pwd)

reg_script="$DIR"/"prepareRegulation.sh"
gwas_script="$DIR"/"prepareGWAS.py"
gtex_script="$DIR"/"prepareGTEx.py"

regdir=${regdir/%\/}
out=${out/%\/}

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating output directories"
mkdir -p "$out/gwas"
mkdir -p "$out/gtex"
mkdir -p "$out/regulation"
mkdir -p "$out/temp"

# ----------------------- DOWNLOADING DATA --------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx data"
wget --quiet -c "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar" -O "$out/temp/GTEx.tar"

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx lookup data"
wget --quiet -c "https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" -O "$out/gtex/GTEx.lookup.txt.gz"

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx Gencode data"
wget --quiet -c "https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf" -O "$out/gtex/gencode.gtf"

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GWAS catalog"
wget --quiet -c "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative" -O "$out/gwas/gwas_full.tsv.gz"

cd "$out/temp"
mkdir -p regulation
cd regulation
echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading Ensembl Regulation"
lftp -c 'open ftp://ftp.ensembl.org/pub/release-100/regulation/homo_sapiens/RegulatoryFeatureActivity/; mirror -P . .' > /dev/null 2> /dev/null

cd "$CDIR"

# -------------------------------------------------------------------------

regdir="$out/temp/regulation"

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating regulation file"
"$reg_script" -i "$regdir" -o "$out/regulation/regulation.bed"

# -------------------------------------------------------------------------

gtex="$out/temp/GTEx.tar"

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GTEx file"
"$gtex_script" -i "$gtex" -g "$out/gtex/gencode.gtf" -l "$out/gtex/GTEx.lookup.txt.gz" | sort -k1,1 -k2,2n > "$out/gtex/gtex.bed"
echo $(date '+%d/%m/%Y %H:%M:%S') "Compressing GTEx file"
bgzip -f "$out/gtex/gtex.bed"
echo $(date '+%d/%m/%Y %H:%M:%S') "Indexing GTEx file"
tabix -f -p bed "$out/gtex/gtex.bed.gz"

# -------------------------------------------------------------------------

gwas="$out/gwas/gwas_full.tsv.gz"

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GWAS file"
PYTHONPATH=$(dirname "$DIR") "$gwas_script" -i "$gwas" | gzip - > "$out/gwas/gwas.tsv.gz"

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Deleting temporary data"
rm -rf "$out/temp/*"

exit 0
