#!/bin/bash


# cat experimentList.tab | awk '($2=="mm10") && ($3=="TFs") {print $0}' > resources/experimentList_mm10_TFs_and_others.tab

# filter metadata
metadata='resources/experimentList_mm10_TFs.tab'
#filter_out=("Epitope tags" "Succinyllysine" "Propionyllysine" "O-GlcNAc" "Lysin homocysteine" "Crotonyl lysine" "Butyryllysine" "8-Hydroxydeoxyguanosine" "5-hmC" "5-mC" "ADP-ribose")
#cat '../resources/experimentList_mm10_TFs_and_others.tab' | \
#    grep -v "Epitope tags" | \
#    grep -v "Succinyllysine" | \
#    grep -v "Propionyllysine" | \
#    grep -v "O-GlcNAc" | \
#    grep -v "Lysin homocysteine" | \
#    grep -v "Crotonyl lysine" | \
#    grep -v "Butyryllysine" | \
#    grep -v "8-Hydroxydeoxyguanosine" | \
#    grep -v "5-hmC" | \
#    grep -v "5-mC" | \
#    grep -v "ADP-ribose" > $metadata


Genome=mm10
Threshold=05
mapfile -t EXPERIMENT_ID  < <(cut -f1 "$metadata")

# get files
for Experimental_ID in "${EXPERIMENT_ID[@]}"
do
    outfile="resources/tracks/${Experimental_ID}.bw"
    if ! [ -e "$outfile" ]
    then
        echo "${Experimental_ID} bigwig"
        
        # BigWig Download URL:
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bw/${Experimental_ID}.bw -O "$outfile"
    fi

    outfile="resources/tracks/${Experimental_ID}.${Threshold}.bb"
    if ! [ -e "$outfile" ]
    then
        echo "${Experimental_ID} bigbed"

        # Peak-call (BigBed) Download URL: #(Threshold = 05, 10, or 20)
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bb${Threshold}/${Experimental_ID}.${Threshold}.bb -O "$outfile"
    fi
done
