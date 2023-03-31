#!/bin/bash

data_file=resources/complex/MusMusculus.tsv
mapfile -t ID < <(cut -f1 $data_file| tail -n+2)
for id in "${ID[@]}"
do
    wget https://www.ebi.ac.uk/intact/complex-ws/export/"$id" -O "resources/complex/json/${id}.json"
done
