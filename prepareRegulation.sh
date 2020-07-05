#!/bin/bash

#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -i <input dir>"
    echo "          -o <output file name>"
    exit 0
}

OPTIND=1
while getopts "i:o:" optname; do
    case "$optname" in
        "i" ) indir="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

for f in $(find "$indir" -maxdepth 1 -name "*.gff.gz");do 
    b=$(basename $f)
    cell_type=$(echo $b| perl -lne '$x="NA";if (/homo_sapiens\.GRCh38\.(\w+)\.Regulatory_Build\.regulatory_activity\.\d+\.gff\.gz/){$x=\1;} print $x;')
    if [[ "$cell_type" != "NA" ]];then
	zcat "$f" | perl -F"\t" -lane '' >> "$out"
    else	
	echo "WARNING: file $f has wrong name pattern"
	echo "Expected pattern: homo_sapiens.GRCh38.\"cell_type_name\".Regulatory_Build.regulatory_activity.\\d+.gff.gz"
	echo "Skipping"
	echo ""
	continue
    fi

