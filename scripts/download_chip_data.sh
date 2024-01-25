#!/bin/bash

# filter metadata
# mm10
# cat experimentList.tab | awk '($2=="mm10") && ($3=="TFs") {print $0}' | \
#     grep -v "Epitope tags" | \
#     grep -v "Succinyllysine" | \
#     grep -v "Propionyllysine" | \
#     grep -v "O-GlcNAc" | \
#     grep -v "Lysin homocysteine" | \
#     grep -v "Crotonyl lysine" | \
#     grep -v "Butyryllysine" | \
#     grep -v "8-Hydroxydeoxyguanosine" | \
#     grep -v "5-hmC" | \
#     grep -v "5-mC" | \
#     grep -v "ADP-ribose" > experimentList_mm10_TFs.tab
# hg38
# cat experimentList.tab | awk '($2=="hg38") && ($3=="TFs") {print $0}' | \
#     grep -v "Pan-acetyllysine" | \
#     grep -v "O-GlcNAc" | \
#     grep -v "MethylCap" | \
#     grep -v "Hepatitis B Virus X antigen" | \
#     grep -v "Epitope tags" | \
#     grep -v "Cyclobutane pyrimidine dimers" | \
#     grep -v "Crotonyllysine" | \
#     grep -v "8-Hydroxydeoxyguanosine" | \
#     grep -v "5-hmC" | \
#     grep -v "5-mC" > experimentList_hg38_TFs.tab
# hg19
# cat experimentList.tab | awk '($2=="hg19") && ($3=="TFs") {print $0}' | \
#     grep -v "Pan-acetyllysine" | \
#     grep -v "O-GlcNAc" | \
#     grep -v "MethylCap" | \
#     grep -v "Hepatitis B Virus X antigen" | \
#     grep -v "Epitope tags" | \
#     grep -v "Cyclobutane pyrimidine dimers" | \
#     grep -v "Crotonyllysine" | \
#     grep -v "8-Hydroxydeoxyguanosine" | \
#     grep -v "5-hmC" | \
#     grep -v "5-mC" > experimentList_hg19_TFs.tab

Genome=hg38
metadata="resources/experimentList_${Genome}_TFs_only_QC_filtered.tab"

mapfile -t EXPERIMENT_ID  < <(cut -f1 "$metadata")
Threshold=05

# get files
for Experimental_ID in "${EXPERIMENT_ID[@]}"
do
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
