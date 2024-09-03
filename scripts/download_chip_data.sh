#!/bin/bash

Genome=$1
metadata=$2

mapfile -t EXPERIMENT_ID  < <(cut -f1 "$metadata")
Threshold=05

# get files
for Experimental_ID in "${EXPERIMENT_ID[@]}"
do
    echo 
    outfile="resources/tracks/${Genome}/${Experimental_ID}.bw"
    if ! [ -e "$outfile" ]
    then
        echo "${Experimental_ID} bigwig"
        
        # BigWig Download URL:
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bw/${Experimental_ID}.bw -O "$outfile"
    fi

    outfile="resources/tracks/${Genome}/${Experimental_ID}.${Threshold}.bb"
    if ! [ -e "$outfile" ]
    then
        echo "${Experimental_ID} bigbed"

        # Peak-call (BigBed) Download URL: #(Threshold = 05, 10, or 20)
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bb${Threshold}/${Experimental_ID}.${Threshold}.bb -O "$outfile"
    fi
done
# Or rsync between servers