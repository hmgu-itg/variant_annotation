#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -o <output dir>"
    exit 0
}

OPTIND=1
while getopts "o:" optname; do
    case "$optname" in
        "o" ) out="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

out=${out/%\/}

mkdir -p $out

tmpfile1=$(mktemp -p $out experiments.XXXX)
tmpfile2=$(mktemp -p $out fnames.XXXX)

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading list of baseline experiments"
getBaselineExperiments.py > $tmpfile1

if [[ $? != 0 ]];then
    echo $(date '+%d/%m/%Y %H:%M:%S') "Error occured while executing getBaselineExperiments.py"
    rm -f $tmpfile1
    rm -f $tmpfile2
    exit 1
fi

echo $(date '+%d/%m/%Y %H:%M:%S') "Downloading data"
cat $tmpfile1 | sed 's/ class="[^"]*"//g' | sed 's/<\/div><\/div>/\n/g'|sed 's/<div><div><span>/\n<div><div><span>/g' |perl -lne 'print $1 if /^\<div\>\<div\>\<span\>(.*)\<\/div\>\<\/label\>$/;' |perl -lne 'print $1."\t".$2."\t".$3 if /.*\/experiments\/(.*)?\/Results.*Design\"\>\<ul\>(.*)?\<\/ul\>\<\/a\>.*\<div\>\<span\>\<ul\>(.*)?\<\/ul\>\<\/span\>\<\/div\>/;' | sed 's/<\/li><li>/;/g' | sed 's/<li>//g' | sed 's/<\/li>//g' | grep "organism part" | grep -v disease | grep -v "cell line" | grep -v individual | cut -f 1,3| tr ' ' '_'| while read id p;do if [[ $p =~ RNA ]];then wget -q https://www.ebi.ac.uk/gxa/experiments-content/"$id"/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv -O "$out"/"$id".tsv;elif [[ $p =~ Proteomics ]];then wget -q https://www.ebi.ac.uk/gxa/experiments-content/"$id"/resources/ExperimentDownloadSupplier.Proteomics/tsv -O "$out"/"$id".tsv;fi;fname="$out"/"$id".tsv;echo $fname >> $tmpfile2;done

echo $(date '+%d/%m/%Y %H:%M:%S') "Concatenating"
concat.py $(cat $tmpfile2| tr '\n' ' ') | gzip - > $out/GXA.tsv.gz

echo $(date '+%d/%m/%Y %H:%M:%S') "Removing files"
cat $tmpfile2 | while read f;do rm -f $f;done

rm -f $tmpfile1
rm -f $tmpfile2

exit 0

