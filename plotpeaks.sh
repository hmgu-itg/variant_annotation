#!/usr/bin/env bash

###################################################################
#
#                          WORKFLOW:
#
# 1) use selectPeaks.py to select peaks from association file
# 2) loop over all peaks and create pdf/html files for each peak
# 2.1)  get peak's chr/pos/id2/id2; for a peak at 1:12345 with alleles a1/a2 = A/AG, id1=1_12345_A_AG, id2=1_12345_AG_A
# 2.2)  for each input PLINK dataset select variants from the <flank_bp> neighbourhood; merge all selected datasets
# 2.3)  select neighbouring variants from association file and save them in "peakdata"
# 2.4)  determine a set of variant common to peakdata and merged
# 2.5)  keep only common variants in peakdata and merged
# 2.6)  use annotateIDList.py to create table with rs IDs, VEP and phenotype annotations
# 2.7)  assign variants in peakdata and merged rsIDs, if possible
# 2.8)  calculate LD
# 2.9)  create locuszoom DB and call locuszoom to create PDF
# 2.10) call interactive_manh.py to create HTML
#
###################################################################


function getColNum () {
    local cmd=$1
    local fname=$2
    local colname=$3
    echo $(fgrep -w $colname <(paste <(seq 1 $($cmd $fname | head -n 1 | tr '\t' '\n'| wc -l)) <($cmd $fname | head -n 1 | tr '\t' '\n')) | cut -f 1)
}

function usage () {
  echo "Usage: plotpeaks -t <p-value threshold>"
  echo "                 -i <input association file>"
  echo "                 -c <chromosome column name>"
  echo "                 -p <p-value column name>"
  echo "                 -m <MAF column name>"
  echo "                 --id <variant ID column name>"
  echo "                 --a1 <allele 1 column name>"
  echo "                 --a2 <allele 2 column name>"
  echo "                 --pos <variant position column name>"
  echo "                 --plink <comma separated PLINK prefixes>"
  echo "                 --flank <flank size (bp)>: optional, default: 500000"
  echo "                 --dbsnp <dbSNP VCF>: optional"
}

signif=""
chrcol=""
pscol=""
rscol=""
mafcol=""
rscol=""
a1col=""
a2col=""
files=""

# default
dbsnp=""
flank_bp=500000

OPTS=$(getopt -o t:i:c:p:m: --long id:a1:a2:pos:plink:flank:dbsnp: -n 'plotpeaks')

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; usage ; exit 1 ; fi

echo "OPTS : $OPTS"
eval set -- "$OPTS"

while true; do
  case "$1" in
    -t ) signif=$2; shift 2 ;;
    -i ) assocfile=$2; shift 2 ;;
    -c ) chrcol=$2; shift 2 ;;
    -p ) pvalcol=$2; shift 2 ;;
    -m ) mafcol=$2; shift 2 ;;
    --id ) rscol=$2; shift 2 ;;
    --a1 ) a1col=$2; shift 2 ;;
    --a2 ) a2col=$2; shift 2 ;;
    --pos ) pscol=$2; shift 2 ;;
    --plink ) files=$2; shift 2 ;;
    --flank ) flank_bp=$2; shift 2 ;;
    --dbsnp ) dbsnp=$2; shift 2 ;;
    * ) echo "Unexpected option: $1" ; usage ; break ;;
  esac
done

echo ""
echo "==================================================="
echo "P-VALUE THRESHOLD: $signif"
echo "ASSOCIATION FILE: $assocfile"
echo "CHROM COLUMN NAME: $chrcol"
echo "POS COLUMN NAME: $pscol"
echo "PVAL COLUMN NAME: $pvalcol"
echo "ID COLUMN NAME: $rscol"
echo "ALLELE 1 COLUMN NAME: $a1col"
echo "ALLELE 2 COLUMN NAME: $a2col"
echo "MAF COLUMN NAME: $mafcol"
echo "FLANK SIZE: $flank_bp"
echo "PLINK FILES: $files"
echo "DBSNP VCF: $dbsnp"
echo "==================================================="
echo ""

if [[ -z "$signif" || -z "$assocfile" || -z "$chrcol" || -z "$pvalcol" || -z "$mafcol" || -z "$rscol" || -z "$a1col" || -z "$a2col" || -z "$pscol" || -z "$files" ]];then
    echo "ERROR: some input parameters were not specified"
    usage
    exit 1
fi

declare -a files=($(echo $files | tr ',' ' '))
declare -a ids=($(seq 0 $(($(echo $files | tr ',' '\n'| wc -l)-1)) | tr '\n' ' '))

flank_kb=$(echo "$flank_bp/1000" | bc)
ext_flank_kb=$((flank_kb+100))
if [[ "$assocfile"~/gz$/ ]];then
	cat=zcat
else
	cat=cat
fi

# get column numbers for all necessary columns
chrcoli=$(getColNum $cat $assocfile $chrcol)
pscoli=$(getColNum $cat $assocfile $pscol)
rscoli=$(getColNum $cat $assocfile $rscol)
pvalcoli=$(getColNum $cat $assocfile $pvalcol)
a1coli=$(getColNum $cat $assocfile $a1col)
a2coli=$(getColNum $cat $assocfile $a2col)

tmp_outdir=$(mktemp -d -p $(pwd) temp_plotpeaks_XXXXXXXX)
if [[ -z "$tmp_outdir" ]];then
    echo "ERROR: failed to create temporary output dir"
    exit 1
fi

echo -e "--------------------------- SELECTING PEAKS -------------------------------\n"
selectPeaks.py -i "$assocfile" -c "$chrcol" -p "$pscol" -v "$pvalcol" -f "$flank_bp" -t "$signif" -o "$tmp_outdir"
if [[ $? -ne 0 ]];then
    echo "ERROR: selectPeaks.py returned non-zero status"
    exit 1
fi
echo ""
echo -e "----------------------------------------------------------------------------\n"

for fname in $(find "$tmp_outdir" -name "peaks_*.txt" | sort);do
    echo -e "INFO: current peak file: $fname\n"
    tail -n +2 "$fname" | while read peakline;do
	echo -e "INFO: current peak line: $peakline\n" 
	read -r peak_chr peak_pos refvar1 refvar2 <<<$(echo "$peakline" | awk -v i="$chrcoli" -v j="$pscoli" -v a="$a1coli" -v b="$a2coli" '{x=toupper($a);y=toupper($b);print $i,$j,$i"_"$j"_"x"_"y,$i"_"$j"_"y"_"x;}')

	echo ""
	echo "peak chr: $peak_chr"
	echo "peak pos: $peak_pos"
	echo "peak ID1: $refvar1"
	echo "peak ID2: $refvar2"
	echo ""

	start_bp=$((peak_pos-flank_bp))
	if [[ $start_bp -lt 1 ]];then
	    start_bp=1
	fi
	end_bp=$((peak_pos+flank_bp))

	# PLINK selecting and merging
	# remove possible "chr" prefix from chromosome names
	# set proper variant IDs in .bim
	rm -f "$tmp_outdir"/mergelist
	for id in "${ids[@]}"
	do
	    echo "Selecting variants from ${files[$id]} using PLINK"
	    plink --memory 15000 --bfile ${files[$id]} --chr $peak_chr --from-bp $start_bp --to-bp $end_bp --out "$tmp_outdir"/$id --make-bed
	    sed 's/^chr//' "$tmp_outdir"/$id.bim | awk 'BEGIN{FS="\t";OFS="\t";}{$2=$1"_"$4"_"$5"_"$6; print $0;}' | sponge "$tmp_outdir"/$id.bim
	    echo "$tmp_outdir"/"$id" >> "$tmp_outdir"/mergelist
	done
	echo -e "PLINK done\n"
	
	echo "Merging"
	plink --merge-list "$tmp_outdir"/mergelist --make-bed --out "$tmp_outdir"/"$peak_chr"."$peak_pos".merged --allow-no-sex
	if [[ -f "$tmp_outdir"/"$peak_chr"."$peak_pos".merged-merge.missnp ]];then
	    for id in "${ids[@]}"
	    do 
		plink --bfile "$tmp_outdir"/$id --exclude "$tmp_outdir"/"$peak_chr"."$peak_pos".merged-merge.missnp --make-bed --out "$tmp_outdir"/$id.tmp
		mv "$tmp_outdir"/$id.tmp.bed "$tmp_outdir"/$id.bed
		mv "$tmp_outdir"/$id.tmp.bim "$tmp_outdir"/$id.bim
		mv "$tmp_outdir"/$id.tmp.fam "$tmp_outdir"/$id.fam
	    done
	    plink --merge-list "$tmp_outdir"/mergelist --make-bed --out "$tmp_outdir"/"$peak_chr"."$peak_pos".merged --allow-no-sex
	fi
	
	# setting refsnp
	if grep -q "$refvar1" "$tmp_outdir"/"$peak_chr"."$peak_pos".merged.bim;then
	    refsnp="$refvar1"
	else
	    if grep -q "$refvar2" "$tmp_outdir"/"$peak_chr"."$peak_pos".merged.bim;then
		refsnp="$refvar2"
	    else
		echo "WARNING: $refvar1 / $refvar2 not found in $tmp_outdir/$peak_chr.$peak_pos.merged.bim; skipping"
		continue
	    fi
	fi
	echo "refsnp=$refsnp"
	echo -e "Done\n"
	
	# select neighbouring variants from association file
	# no header in peakdata
	echo "Selecting neighbouring variants from $assocfile: ${peak_chr}:${start_bp}-${end_bp}"
	$cat $assocfile | awk -v i="$chrcoli" -v c="$peak_chr" -v j="$pscoli" -v p="$peak_pos" -v f="$flank_bp" 'BEGIN{FS="\t";OFS="\t";}{if ($i==c && $j>(p-f) && $j<(p+f)){print $0;}}' | grep -v nan > "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata
	echo -e "Done\n"

	# common variants between merged PLINK and peakdata
	echo "Extracting common variants between peakdata and merged file"
	cat <(cat "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | awk -v i="$chrcoli" -v j="$pscoli" -v a="$a1coli" -v b="$a2coli" '{x=toupper($a);y=toupper($b);print $i"_"$j"_"x"_"y; print $i"_"$j"_"y"_"x;}') <(cut -f 2 "$tmp_outdir"/"$peak_chr"."$peak_pos".merged.bim) | sort| uniq -d > "$tmp_outdir"/"$peak_chr"."$peak_pos".common
	echo -e "Done\n"
	
	# setting variant IDs in peakdata to those in "common" file
	# after setting IDs, peakdata may contain less variants than "merged_common" PLINK
	# because "merged" PLINK and "common" may contain variants
	# with same position and permuted alleles: 1_12345_A_C and 1_12345_C_A
	# whereas peakdata contains only one of those (the first match in "common")
	# this means LD data (based on PLINK data) may contain more variants than peakdata.chrpos
	# this later leads to locuszoom complaining:
	# Warning: could not find position for SNP 22_49623908_G_A in user-supplied --ld file, skipping..
	
	echo "Setting variant IDs in" "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata
	rm -f "$tmp_outdir"/missing
	cat "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | while read line;do read -r var1 var2 <<<$(echo "$line" | awk -v i="$chrcoli" -v j="$pscoli" -v a="$a1coli" -v b="$a2coli" '{x=toupper($a);y=toupper($b);print $i"_"$j"_"x"_"y,$i"_"$j"_"y"_"x;}'); if grep -w -q "$var1" "$tmp_outdir"/"$peak_chr"."$peak_pos".common; then ID="$var1";else if grep -w -q "$var2" "$tmp_outdir"/"$peak_chr"."$peak_pos".common; then ID="$var2";else echo "$line" >>"$tmp_outdir"/"$peak_chr"."$peak_pos".missing;continue;fi;fi;echo "$line" | awk -v i="$rscoli" -v x="$ID" '{$i=x;print $0;}';done | tr ' ' '\t' | sponge "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata
	echo -e "Done\n"

	# selecting common variants from merged
	echo "Selecting common variants from $peak_chr.$peak_pos.merged"
	plink --bfile "$tmp_outdir"/"$peak_chr"."$peak_pos".merged --extract "$tmp_outdir"/"$peak_chr"."$peak_pos".common --allow-no-sex --make-bed --out "$tmp_outdir"/mergedtemp
	mv "$tmp_outdir"/mergedtemp.bed "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_common.bed
	mv "$tmp_outdir"/mergedtemp.bim "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_common.bim
	mv "$tmp_outdir"/mergedtemp.fam "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_common.fam
	echo -e "Done\n"
	
	# create table with rs IDs, VEP and phenotype annotations
	echo "Annotating peakdata variants"
	if [[ -z "$dbsnp" ]];then
	    cut -f $rscoli  "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | annotateIDList.py -v debug 1>"$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table 2>"$tmp_outdir"/"$peak_chr"."$peak_pos".debug
	else
	    minp=$(cut -f 2 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | sort -n | head -n 1)
	    maxp=$(cut -f 2 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | sort -nr | head -n 1)
	    tabix "$dbsnp" "$peak_chr":"$minp"-"$maxp" | cut -f 2-5 > "$tmp_outdir"/"$peak_chr"."$peak_pos".dbsnp
	    join -a 2 -j 1 "$tmp_outdir"/"$peak_chr"."$peak_pos".dbsnp <(cut -f 2,3 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | sort -k1,1n) | awk 'NF==2{print $2;}' > "$tmp_outdir"/"$peak_chr"."$peak_pos".not_in_dbsnp
	    # 22 16059596 rs1317973428 GA G is the same as 22_16059596_GAA_GAAA
	    join -a 2 -j 1 "$tmp_outdir"/"$peak_chr"."$peak_pos".dbsnp <(cut -f 2,3 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | sort -k1,1n) | awk 'NF!=2'| compareVariants.pl 1>"$tmp_outdir"/"$peak_chr"."$peak_pos".in_dbsnp 2>>"$tmp_outdir"/"$peak_chr"."$peak_pos".not_in_dbsnp
	    cat "$tmp_outdir"/"$peak_chr"."$peak_pos".not_in_dbsnp | annotateIDList.py -v debug 1>"$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table 2>"$tmp_outdir"/"$peak_chr"."$peak_pos".debug
	    cat "$tmp_outdir"/"$peak_chr"."$peak_pos".in_dbsnp | annotateIDList.py --rs -v debug 1>>"$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table 2>>"$tmp_outdir"/"$peak_chr"."$peak_pos".debug
	fi
	echo -e "Done\n"

	# rename variants in peakdata and merged.bim
	echo "Renaming variants using annotated table"
	join -1 1 -2 $rscoli <(cut -f 1,2 "$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table | sort -k1,1) <(sort -k"$rscoli","$rscoli" "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata) | cut -d ' ' -f 2- | awk 'BEGIN{OFS="\t";}{x=$1;$1=$2;$2=$3;$3=x;print $0;}' | sponge "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata
	cat  <(echo -e "snp\tchr\tpos") <(cat "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | awk -v i=$rscoli -v j=$chrcoli -v k=$pscoli 'BEGIN{FS="\t";OFS="\t";}{print $i,$j,$k;}') > "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata.chrpos
	cut -f 1,2 "$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table > "$tmp_outdir"/"$peak_chr"."$peak_pos".rename
	plink --bfile "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_common --out "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_rn --make-bed --allow-no-sex --update-name "$tmp_outdir"/"$peak_chr"."$peak_pos".rename
	echo -e "Done\n"

	# add header to peakdata
	cat <($cat $assocfile| head -n 1) <(cat "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata) | sponge "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata
	
	# update refsnp ID
	refsnp=$(grep -w $refsnp "$tmp_outdir"/"$peak_chr"."$peak_pos".rename| cut -f 2)
	echo "updated refsnp=$refsnp"
	
	# compute LD
	echo "Computing LD"
	plink --bfile "$tmp_outdir"/"$peak_chr"."$peak_pos".merged_rn --r2 dprime --ld-snp $refsnp --ld-window-kb $ext_flank_kb --ld-window 999999 --ld-window-r2 0 --out "$tmp_outdir"/merged
	cat <(echo "snp1 snp2 dprime rsquare") <(tail -n +2 "$tmp_outdir"/merged.ld| sed -e 's/^  *//' -e 's/  */ /g'| awk '{print $6,$3,$7,$8}') > "$tmp_outdir"/"$peak_chr"."$peak_pos".ld
	echo -e "Done\n"

	# create locuszoom DB for the current peak neighborhood
	echo "Creating locuszoom DB"
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --snp_pos "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata.chrpos
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --refflat /opt/locuszoom/data/refFlat_b38.txt
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --recomb_rate /opt/locuszoom/data/recomb_rate_b38.txt
	echo -e "Done\n"

	# call locuszoom
	echo "Calling locuszoom"
	locuszoom --metal "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db "$tmp_outdir"/locuszoom.db --prefix ${peak_chr}.${peak_pos}.${flank_bp} --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld "$tmp_outdir"/"$peak_chr"."$peak_pos".ld --start="$start_bp" --end="$end_bp" --chr="$peak_chr" showRecomb=T --build b38
	echo -e "Done\n"

	# reformat traits column
	cat "$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table | perl -lne '@a=split(/\t/);@b=split(/\s*:\s*/,$a[2]);$b[0]=~s/[{"]//g;$b[1]=~s/[}"]//g;$a[3]=~s/[][]//g;$a[3]="NA" if $a[3]=~/^\s*$/;$,="\t";$a[3]=~s/,\s+/,/g;$a[3]=~s/\s+/_/g;print $a[0],$a[1],$b[0],$b[1],$a[3];' > "$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table_mod 
	
	# prepare data for interactive manhattan plotting
	echo "Joining"
	join --header -1 $rscoli -2 1 <(cat <(head -n 1 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata) <(tail -n +2 "$tmp_outdir"/"$peak_chr"."$peak_pos".peakdata | sort -k$rscoli,$rscoli)) <(cat <(echo $rscol ld) <(tail -n +2 "$tmp_outdir"/"$peak_chr"."$peak_pos".ld | tr ' ' '\t' | cut -f 1,3 | sort -k1,1)) | tr ' ' '\t' > "$tmp_outdir"/"$peak_chr"."$peak_pos".join1
	join --header -1 1 -2 1 <(cat <(head -n 1 "$tmp_outdir"/"$peak_chr"."$peak_pos".join1) <(sort -k1,1 "$tmp_outdir"/"$peak_chr"."$peak_pos".join1)) <(cat <(echo $rscol gene consequence traits) <(cut -f 2- "$tmp_outdir"/"$peak_chr"."$peak_pos".annotated_table_mod | sort -k1,1)) | tr ' ' '\t' | sed -e 's/"//g' -e 's/_/ /g' > "$tmp_outdir"/"$peak_chr"."$peak_pos".join2
	echo -e "Done\n"

	# create ineractive HTML
	echo "Creating HTML"
	#PYTHONPATH=~/variant_annotation/python/ ~/variant_annotation/interactive_manh.py "$tmp_outdir"/"$peak_chr"."$peak_pos".join2 $chrcol $pvalcol $pscol $rscol $mafcol "$peak_chr"."$peak_pos"."$flank_bp".html # <-- change this line
	interactive_manh.py "$tmp_outdir"/"$peak_chr"."$peak_pos".join2 $chrcol $pvalcol $pscol $rscol $mafcol "$peak_chr"."$peak_pos"."$flank_bp".html
	echo -e "Done\n"
    done
done

rm -rf "$tmp_outdir"

exit 0

