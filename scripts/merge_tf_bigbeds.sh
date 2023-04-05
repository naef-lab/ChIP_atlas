#!/bin/bash

mapfile -t TF < <(cut -f4 resources/experimentList_mm10_TFs.tab | sort -u)

for tf in "${TF[@]}"
do
    echo "$tf"

    # define output and tmp bedfile for tf and empty if exists
    fout=results/TF_bedfiles/${tf}.bed
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
    mapfile -t ID < <(awk -F '\t' -v tf="$tf" '$4 == tf {print $1}' resources/experimentList_mm10_TFs.tab)
    for id in "${ID[@]}"
    do
        # get input file
        fin=resources/tracks/${id}.05.bb

        # append bed to outfile
        bigBedToBed "$fin" tmp/tmp.bed
        cat tmp/tmp.bed >> "$tmp_fout"
    done
    # sort and merge
    sort -k1,1 -k2,2n "$tmp_fout" | bedtools merge -c 4 -o max > "$fout"

    # TODO to bigbeds
    # fout_bb=results/bedfiles/${tf}.bb
    # bedToBigBed in.bed chrom.sizes out.bb
done
rm tmp/tmp.bed