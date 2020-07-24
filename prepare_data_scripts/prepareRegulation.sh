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

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

#tmpfile=$(mktemp "$PWD"/tmp_regulation.XXXX)

for f in $(find "$indir" -maxdepth 1 -name "*.gff.gz");do 
    b=$(basename $f)
    #echo "Current file: $b"
    cell_type=$(echo $b| perl -lne '$x="NA";if (/homo_sapiens\.GRCh38\.(\w+)\.Regulatory_Build\.regulatory_activity\.\d+\.gff\.gz/){$x=$1;} print $x;')
    #echo "Cell type: $cell_type"
    if [[ "$cell_type" != "NA" ]];then
	export cell_type="$cell_type"
	zcat "$f" | perl -F"\t" -lane '$ct=$ENV{cell_type};$,="\t";$a="NA";$id="NA";if ($F[8] =~ /activity=([^;]+);/){$a=$1;}if ($F[8] =~ /regulatory_feature_stable_id=([^;]+)/){$id=$1;} print $F[0],$F[3]-1,$F[4],$F[2],$ct,$a,$id;' >> "$out"
    else	
	echo "WARNING: file $f has wrong name pattern"
	echo "Expected pattern: homo_sapiens.GRCh38.\"cell_type_name\".Regulatory_Build.regulatory_activity.\\d+.gff.gz"
	echo "Skipping"
	echo ""
	continue
    fi
done

echo "Sorting"
sort -k1,1 -k2,2n "$out" | sponge "$out"
echo "Compressing"
bgzip -f "$out"
echo "Indexing"
tabix -p bed "$out".gz
#rm "$tmpfile"
echo "Done"

