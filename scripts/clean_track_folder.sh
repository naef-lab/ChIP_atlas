#!/bin/bash

track_folder="resources/tracks"
GENOME=(hg38 mm10)

for genome in "${GENOME[@]}"
do
    experiment_list="resources/experimentList_${genome}_TFs_only_QC_filtered.tab"
    echo "$genome"

    n_tot=$(find "$track_folder"/"$genome"/*.bw|wc -l)
    n_in=0
    n_out=0
    for id in $(find "$track_folder"/"$genome"/*.bw|cut -f4 -d'/'|cut -f1 -d'.')
    do
        n=$(grep -c "^$id\>" "$experiment_list")
        if [ "$n" -eq 0 ]
        then
            n_out=$((n_out+1))
            rm "$track_folder"/"$genome"/"$id".*
        else
            n_in=$((n_in+1))
        fi
    done

    echo "n_tot: $n_tot"
    echo "n_in: $n_in"
    echo "n_exp: $(wc -l "$experiment_list")"
    echo "n_out: $n_out"
done