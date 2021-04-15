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

signif=$1
assocfile=$2
chrcol=$3
pscol=$4
rscol=$5
pvalcol=$6
a1col=$7
a2col=$8
mafcol=$9
files=${10}
flank_bp=${11}
filelist=$files

echo
echo
echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}
echo
echo

declare -a files=($(echo $files | tr ',' ' '))
declare -a ids=($(seq 0 $(($(echo ${10} | tr ',' '\n'| wc -l)-1)) | tr '\n' ' '))
if [ -z "$flank_bp" ]; then
	flank_bp=500000
fi

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

echo COLUMNS "CHR: "$chrcoli . "POS: "$pscoli . "rsID: "$rscoli . "Pval: "$pvalcoli . "A1: "$a1coli . "A2: "$a2coli

echo "Looking for peaks..."
echo "===================="
echo "(p-value $pvalcol - $pvalcoli; $signif ; $assocfile ; $cat)"
echo
echo

tmp_outdir=$(mktemp -d -p $(pwd) temp_plotpeaks_XXXXXXXX)
if [[ -z "$tmp_outdir" ]];then
    echo "ERROR: failed to create temporary output dir"
    exit 1
fi

echo -e "--------------------------- SELECTING PEAKS -------------------------------\n"
~/variant_annotation/selectPeaks.py -i "$assocfile" -c "$chrcol" -p "$pscol" -v "$pvalcol" -f "$flank_bp" -t "$signif" -o "$tmp_outdir" # <--- change this line
if [[ $? -ne 0 ]];then
    echo "ERROR: selectPeaks.py returned non-zero status"
    exit 1
fi
echo ""
echo -e "----------------------------------------------------------------------------\n"

for fname in $(find "$tmp_outdir" -name "peak*chr21*.txt" | sort);do # <--- change this line
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
	# also remove possible "chr" prefix from variant IDs and chromosome names
	rm -f "$tmp_outdir"/mergelist
	for id in "${ids[@]}"
	do
	    echo "Selecting variants from ${files[$id]} using PLINK"
	    plink --memory 15000 --bfile ${files[$id]} --chr $peak_chr --from-bp $start_bp --to-bp $end_bp --out "$tmp_outdir"/$id --make-bed
	    cat "$tmp_outdir"/$id.bim | sed 's/chr//g' | sponge "$tmp_outdir"/$id.bim
	    echo "$tmp_outdir"/"$id" >> "$tmp_outdir"/mergelist
	done
	echo -e "PLINK done\n"
	
	echo "Merging"
	plink --merge-list "$tmp_outdir"/mergelist --make-bed --out "$tmp_outdir"/merged --allow-no-sex
	if [[ -f "$tmp_outdir"/merged-merge.missnp ]];then
	    #grep -v -w -f "$tmp_outdir"/merged-merge.missnp "$tmp_outdir"/peakdata | sponge "$tmp_outdir"/peakdata 
	    for id in "${ids[@]}"
	    do 
		plink --bfile "$tmp_outdir"/$id --exclude "$tmp_outdir"/merged-merge.missnp --make-bed --out "$tmp_outdir"/$id.tmp
		mv "$tmp_outdir"/$id.tmp.bed "$tmp_outdir"/$id.bed
		mv "$tmp_outdir"/$id.tmp.bim "$tmp_outdir"/$id.bim
		mv "$tmp_outdir"/$id.tmp.fam "$tmp_outdir"/$id.fam
	    done
	    plink --merge-list "$tmp_outdir"/mergelist --make-bed --out "$tmp_outdir"/merged --allow-no-sex
	fi
	
	# setting refsnp
	if grep -q "$refvar1" "$tmp_outdir"/merged.bim;then
	    refsnp="$refvar1"
	else
	    if grep -q "$refvar2" "$tmp_outdir"/merged.bim;then
		refsnp="$refvar2"
	    else
		echo "WARNING: $refvar1 / $refvar2 not found in $tmp_outdir/merged.bim; skipping"
		continue
	    fi
	fi	
	echo -e "Done\n"
	
	# select neighbouring variants from association file
	# no header in peakdata
	echo "Selecting neighbouring variants from $assocfile: ${peak_chr}:${start_bp}-${end_bp}"
	$cat $assocfile | awk -v i="$chrcoli" -v c="$peak_chr" -v j="$pscoli" -v p="$peak_pos" -v f="$flank_bp" 'BEGIN{FS="\t";OFS="\t";}{if ($i==c && $j>(p-f) && $j<(p+f)){print $0;}}' | grep -v nan > "$tmp_outdir"/peakdata # <--- TODO: better filter for NaN p-values ?
	echo -e "Done\n"

	# common variants between merged PLINK and peakdata
	echo "Extracting common variants between peakdata and merged file"
	cat <(cat "$tmp_outdir"/peakdata | awk -v i="$chrcoli" -v j="$pscoli" -v a="$a1coli" -v b="$a2coli" '{x=toupper($a);y=toupper($b);print $i"_"$j"_"x"_"y; print $i"_"$j"_"y"_"x;}') <(cut -f 2 "$tmp_outdir"/merged.bim) | sort| uniq -d > "$tmp_outdir"/common
	echo -e "Done\n"
	
	# selecting common variants from peakdata
	# setting variant IDs in peakdata
	echo "Setting variant IDs"
	rm -f "$tmp_outdir"/missing
	cat "$tmp_outdir"/peakdata | while read line;do read -r var1 var2 <<<$(echo "$line" | awk -v i="$chrcoli" -v j="$pscoli" -v a="$a1coli" -v b="$a2coli" '{x=toupper($a);y=toupper($b);print $i"_"$j"_"x"_"y,$i"_"$j"_"y"_"x;}'); if grep -q -m 1 "$var1" "$tmp_outdir"/common; then ID="$var1";else if grep -q -m 1 "$var2" "$tmp_outdir"/common; then ID="$var2";else echo "$line" >>"$tmp_outdir"/missing;continue;fi;fi;echo "$line" | awk -v i="$rscoli" -v x="$ID" '{$i=x;print $0;}';done | tr ' ' '\t' |  sponge "$tmp_outdir"/peakdata
	echo -e "Done\n"

	# selecting common variants from merged
	echo "Selecting common variants from merged"
	plink --bfile "$tmp_outdir"/merged --extract "$tmp_outdir"/common --allow-no-sex --make-bed --out "$tmp_outdir"/mergedtemp
	mv "$tmp_outdir"/mergedtemp.bed "$tmp_outdir"/merged.bed
	mv "$tmp_outdir"/mergedtemp.bim "$tmp_outdir"/merged.bim
	mv "$tmp_outdir"/mergedtemp.fam "$tmp_outdir"/merged.fam
	echo -e "Done\n"
	
	# create table with rs IDs, VEP and phenotype annotations
	echo "Annotating peakdata IDs"
	cut -f $rscoli  "$tmp_outdir"/peakdata | PYTHONPATH=~/variant_annotation/python/ ~/variant_annotation/annotateIDList.py -v debug 1>"$tmp_outdir"/annotated_table 2>"$tmp_outdir"/debug # <--- change this line
	echo -e "Done\n"

	# rename variants in peakdata and merged.bim
	echo "Renaming variants"
	join -1 1 -2 $rscoli <(cut -f 1,2 "$tmp_outdir"/annotated_table | sort -k1,1) <(sort -k"$rscoli","$rscoli" "$tmp_outdir"/peakdata) | cut -d ' ' -f 2- | awk 'BEGIN{OFS="\t";}{x=$1;$1=$2;$2=$3;$3=x;print $0;}' | sponge "$tmp_outdir"/peakdata
	cat  <(echo -e "snp\tchr\tpos") <(cat "$tmp_outdir"/peakdata | awk -v i=$rscoli -v j=$chrcoli -v k=$pscoli 'BEGIN{FS="\t";OFS="\t";}{print $i,$j,$k;}') > "$tmp_outdir"/peakdata.chrpos
	cut -f 1,2 "$tmp_outdir"/annotated_table > "$tmp_outdir"/rename
	plink --bfile "$tmp_outdir"/merged --out "$tmp_outdir"/mergedtemp --make-bed --allow-no-sex --update-name "$tmp_outdir"/rename
	mv "$tmp_outdir"/mergedtemp.bed "$tmp_outdir"/merged.bed
	mv "$tmp_outdir"/mergedtemp.bim "$tmp_outdir"/merged.bim
	mv "$tmp_outdir"/mergedtemp.fam "$tmp_outdir"/merged.fam	
	echo -e "Done\n"

	# add header to peakdata
	cat <($cat $assocfile| head -n 1) <(cat "$tmp_outdir"/peakdata) | sponge "$tmp_outdir"/peakdata
	
	# update refsnp ID
	refsnp=$(grep -w $refsnp "$tmp_outdir"/rename| cut -f 2)
	
	# compute LD
	echo "Computing LD"
	plink --bfile "$tmp_outdir"/merged --r2 dprime --ld-snp $refsnp --ld-window-kb $ext_flank_kb --ld-window 999999 --ld-window-r2 0 --out "$tmp_outdir"/merged
	cat <(echo "snp1 snp2 dprime rsquare") <(tail -n +2 "$tmp_outdir"/merged.ld| sed -e 's/^  *//' -e 's/  */ /g'| awk '{print $6,$3,$7,$8}') > "$tmp_outdir"/"$peak_chr"."$peak_pos".ld
	echo -e "Done\n"

	# create locuszoom DB for the current peak neighborhood
	echo "Creating locuszoom DB"
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --snp_pos "$tmp_outdir"/peakdata.chrpos
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --refflat /opt/locuszoom/data/refFlat_b38.txt
	dbmeister.py --db "$tmp_outdir"/locuszoom.db --recomb_rate /opt/locuszoom/data/recomb_rate_b38.txt
	echo -e "Done\n"

	# call locuszoom
	echo "Calling locuszoom"
	locuszoom --metal "$tmp_outdir"/peakdata --refsnp "$refsnp" --markercol "$rscol" --pvalcol "$pvalcol" --db "$tmp_outdir"/locuszoom.db --prefix ${peak_chr}.${peak_pos}.${flank_bp} --plotonly showAnnot=T showRefsnpAnnot=T annotPch="21,24,24,25,22,22,8,7" rfrows=20 geneFontSize=.4 --ld "$tmp_outdir"/"$peak_chr"."$peak_pos".ld --start="$start_bp" --end="$end_bp" --chr="$peak_chr" showRecomb=T --build b38
	echo -e "Done\n"

	# modify annotated_table
	cat "$tmp_outdir"/annotated_table | perl -lne '@a=split(/\t/);@b=split(/\s*:\s*/,$a[2]);$b[0]=~s/[{"]//g;$b[1]=~s/[}"]//g;$a[3]=~s/[][]//g;$a[3]="NA" if $a[3]=~/^\s*$/;$,="\t";$a[3]=~s/,\s+/,/g;$a[3]=~s/\s+/_/g;print $a[0],$a[1],$b[0],$b[1],$a[3];' > "$tmp_outdir"/annotated_table_mod 
	
	# prepare data for interactive manhattan plotting
	echo "Joining"
	join --header -1 $rscoli -2 1 <(cat <(head -n 1 "$tmp_outdir"/peakdata) <(tail -n +2 "$tmp_outdir"/peakdata | sort -k$rscoli,$rscoli)) <(cat <(echo $rscol ld) <(tail -n +2 "$tmp_outdir"/"$peak_chr"."$peak_pos".ld | tr ' ' '\t' | cut -f 1,3 | sort -k1,1)) | tr ' ' '\t' > "$tmp_outdir"/join1
	join --header -1 1 -2 1 <(cat <(head -n 1 "$tmp_outdir"/join1) <(sort -k1,1 "$tmp_outdir"/join1)) <(cat <(echo $rscol gene consequence traits) <(cut -f 2- "$tmp_outdir"/annotated_table_mod | sort -k1,1)) | tr ' ' '\t' > "$tmp_outdir"/join2
	echo -e "Done\n"

	# create ineractive HTML

	
    done
done

exit 0

