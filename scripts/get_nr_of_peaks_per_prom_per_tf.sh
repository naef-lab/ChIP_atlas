#!/bin/bash

#promoterome=$1
promoterome=../Promoterome/results/mm10/promoterome_pm2kb_filtered.bed
n_prom=$(wc -l $promoterome | cut -d' ' -f1)

# chip_experiment_list=$2
chip_experiment=resources/experimentList_mm10_TFs_only_QC_filtered.tab
# get unique TFs
mapfile -t TFs < <(cut -f4 $chip_experiment | tail -n+2 | sort -u)

# get genome
# genome=$3
genome='mm10'


# output file
outfile=results/mm10/nr_of_peaks_per_prom_per_tf.tsv

# write header
printf "promoter_tf" > $outfile
for tf in "${TFs[@]}"
do
    printf "\t%s" "$tf" >> $outfile
done
printf "\n" >> $outfile

# loop on promoterome (skip header) and TFs
for n in $(seq 2 "${n_prom}")
do
    echo $n
    id=$(sed "$n""q;d" $promoterome|cut -f6,7|tr '\t' '_')
    printf "%s" "$id" >> $outfile
    for tf in "${TFs[@]}"
    do
        # get nr of peaks
        N=$(bedtools intersect -a results/${genome}/TF_prom_bedfiles/"${tf}".bed -b <(sed "$n"'q;d' $promoterome) -wa | wc -l)
        #bedtools intersect -a results/${genome}/TF_ -b -wa | wc -l
        printf "\t%u" "$N" >> $outfile
    done
    printf "\n" >> $outfile
done
