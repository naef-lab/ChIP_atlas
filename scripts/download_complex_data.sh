#!/bin/bash

# Download complex data from IntAct
# Usage: bash download_complex_data.sh
# Output: resources/complex/json/{complex_id}.json

DATA_FILES=("resources/complex/MusMusculus.tsv" "resources/complex/HomoSapiens.tsv")

for data_file in "${DATA_FILES[@]}"; do

    echo "Downloading data from $data_file"

    mapfile -t ID < <(cut -f1 $data_file| tail -n+2)
    for id in "${ID[@]}"
    do
        if [ ! -f "resources/complex/json/${id}.json" ]; then
            echo "Downloading $id"
            wget https://www.ebi.ac.uk/intact/complex-ws/export/"$id" -O "resources/complex/json/${id}.json"
        fi
    done
done
