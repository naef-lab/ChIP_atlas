#!/bin/bash

genome=$1
tf=$2
chip_experiment=$3
fout=$4

echo "$tf"

# define output and tmp bedfile for tf and empty if exists
if [ -e "$fout" ]
then
    rm "$fout"
fi
touch "$fout"

tmp_fout=tmp/tmp_${tf}.bed
if [ -e "$tmp_fout" ]
then
    rm "$tmp_fout"
fi
touch "$tmp_fout"

# get all experiment ids for tf
mapfile -t ID < <(awk -F '\t' -v tf="$tf" '$4 == tf {print $1}' $chip_experiment)
for id in "${ID[@]}"
do
    # get input file
    fin=resources/tracks/${genome}/${id}.05.bb

    tmp_file=tmp/tmp_${id}.bed

    # append bed to outfile
    bigBedToBed "$fin" "$tmp_file"
    cat "$tmp_file" >> "$tmp_fout"

    rm "$tmp_file"
done
# sort and merge
sort -k1,1 -k2,2n "$tmp_fout" | bedtools merge -c 4 -o max > "$fout"
rm "$tmp_fout"

# TODO to bigbeds
# fout_bb=results/bedfiles/${tf}.bb
# bedToBigBed in.bed chrom.sizes out.bb