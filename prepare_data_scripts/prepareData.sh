#!/bin/bash

function usage {
    echo ""
    echo "Usage:" $(basename $0) "-o <output dir> { -h -g -x -r -w -a}"
    echo " -h : this help message"
    echo " -g -x -r -w -a are optional, for selecting which data to prepare (GWAS, GTEx, Regulation, GWAVA or GXA), by default all."

    exit 0
}

OPTIND=1

prepgwas=0
prepgtex=0
prepreg=0
prepgxa=0
prepgwava=0

while getopts "gxrwao:h" optname; do
    case "$optname" in
        "g" ) prepgwas=1;;
        "x" ) prepgtex=1;;
        "r" ) prepreg=1;;
        "w" ) prepgwava=1;;
        "a" ) prepgxa=1;;
        "o" ) out="${OPTARG}";;
        "h" ) usage ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
fi

if [ "$prepgxa" -eq 0 ] &&[ "$prepgwas" -eq 0 ] && [ "$prepgtex" -eq 0 ] && [ "$prepreg" -eq 0 ] && [ "$prepgwava" -eq 0 ];then
    prepgwas=1
    prepgtex=1
    prepreg=1
    prepgwava=1
    prepgxa=1
fi

pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
  while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`; 
  SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null
DIR=$SCRIPT_PATH
CDIR=$(pwd)

reg_script="$DIR"/"prepareRegulation.sh"
gwas_script="$DIR"/"prepareGWAS.py"
gtex_script="$DIR"/"prepareGTEx.py"
gxa_script="$DIR"/"prepareGXA.sh"

regdir=${regdir/%\/}
out=${out/%\/}

mkdir -p "$out"
logfile="$out"/prepare_data.log
: > "$logfile"

# ------------------------------- URLs ------------------------------------

gtexURL="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar"
gtexlookupURL="https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
gtexgencodeURL="https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"
gwasURL="https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
ensregURL="ftp://ftp.ensembl.org/pub/release-100/regulation/homo_sapiens/RegulatoryFeatureActivity/"
gwavaURL="ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/"
hg19chrom="https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Creating output directories"|tee -a "$logfile"
mkdir -p "$out/gwas"
mkdir -p "$out/gtex"
mkdir -p "$out/regulation"
mkdir -p "$out/temp"
mkdir -p "$out/gwava"
mkdir -p "$out/gxa"

# ----------------------- DOWNLOADING DATA --------------------------------

if [[ "$prepgtex" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx data"|tee -a "$logfile"
    wget --quiet -c "$gtexURL" -O "$out/temp/GTEx.tar"
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $gtexURL"|tee -a "$logfile"
	exit 1
    fi

    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx lookup data"|tee -a "$logfile"
    wget --quiet -c "$gtexlookupURL" -O "$out/gtex/GTEx.lookup.txt.gz"
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $gtexlookupURL"|tee -a "$logfile"
	exit 1
    fi

    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GTEx Gencode data"|tee -a "$logfile"
    wget --quiet -c "$gtexgencodeURL" -O "$out/gtex/gencode.gtf"
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $gtexgencodeURL"|tee -a "$logfile"
	exit 1
    fi
    
    gzip -f "$out/gtex/gencode.gtf"
fi

if [[ "$prepgwas" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GWAS catalog"|tee -a "$logfile"
    wget --quiet -c "$gwasURL" -O "$out/gwas/gwas_full.tsv"
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $gwasURL"|tee -a "$logfile"
	exit 1
    fi
    
    gzip -f "$out/gwas/gwas_full.tsv"
fi

if [[ "$prepreg" -eq 1 ]];then
    cd "$out/temp"
    mkdir -p regulation
    cd regulation
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading Ensembl Regulation"|tee -a "$logfile"
    lftp -c "open $ensregURL; mirror -P . ." >> "$logfile" 2>> "$logfile"
    cd "$CDIR"
fi


if [[ "$prepgwava" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading GWAVA"|tee -a "$logfile"
    wget -r -N -l inf --cut-dirs 5 -np -nH -P $out/gwava --quiet -c "$gwavaURL"
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $gwavaURL"|tee -a "$logfile"
	exit 1
    fi

    wget -q "$hg19chrom" -O $out/gwava/source_data/hg19/human.hg19.genome
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not download $hg19chrom"|tee -a "$logfile"
	exit 1
    fi
fi

# -------------------------- REGULATION -----------------------------------

if [[ "$prepreg" -eq 1 ]];then
    regdir="$out/temp/regulation"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating regulation file"|tee -a "$logfile"
    "$reg_script" -i "$regdir" -o "$out/regulation/regulation.bed" >> "$logfile"
fi

# ----------------------------- GTEx --------------------------------------

if [[ "$prepgtex" -eq 1 ]];then
    gtex="$out/temp/GTEx.tar"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GTEx file"|tee -a "$logfile"
    "$gtex_script" -i "$gtex" -g "$out/gtex/gencode.gtf.gz" -l "$out/gtex/GTEx.lookup.txt.gz" -t "$out/temp" | sort -k1,1 -k2,2n > "$out/gtex/gtex.bed"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Compressing GTEx file"|tee -a "$logfile"
    bgzip -f "$out/gtex/gtex.bed"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Indexing GTEx file"|tee -a "$logfile"
    tabix -f -p bed "$out/gtex/gtex.bed.gz"
    rm "$out/gtex/gencode.gtf.gz" "$out/gtex/GTEx.lookup.txt.gz"
fi

# ----------------------------- GWAS ---------------------------------------

if [[ "$prepgwas" -eq 1 ]];then
    gwas="$out/gwas/gwas_full.tsv.gz"
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GWAS file"|tee -a "$logfile"
    "$gwas_script" -i "$gwas" 2>>"$logfile" | gzip - > "$out/gwas/gwas.tsv.gz" 
fi

# ------------------------ POST-PROCESSING GWAVA --------------------------

if [[ "$prepgwava" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Patching GWAVA"|tee -a "$logfile"
    patch -b $out/gwava/src/gwava_annotate.py /usr/local/bin/variant_annotation/patches/gwava_annotate.py.patch
    patch -b $out/gwava/src/gwava.py /usr/local/bin/variant_annotation/patches/gwava.py.patch
    zcat "$out/gwava/source_data/encode/Gencodev10_TSS_May2012.gff.gz" | sort -k1,1 -k4,5n | gzip - > "$out/temp/tmp.gff.gz"
    mv "$out/temp/tmp.gff.gz" "$out/gwava/source_data/encode/Gencodev10_TSS_May2012.gff.gz"
fi

# ----------------------------- GXA ---------------------------------------

if [[ "$prepgxa" -eq 1 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Creating GXA file"|tee -a "$logfile"
    "$gxa_script" -o "$out/gxa"
fi

# -------------------------------------------------------------------------

echo $(date '+%d/%m/%Y %H:%M:%S') "Deleting temporary data"|tee -a "$logfile"
rm -rf "$out/temp/"

exit 0
