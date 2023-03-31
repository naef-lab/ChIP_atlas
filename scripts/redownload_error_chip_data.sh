#!/bin/bash

infile='resources/to_dowload.txt'
mapfile -t to_download  < <(cat "$infile")

Genome="mm10"
# get files
for Experimental_ID in "${to_download[@]}"
do
    
    outfile="resources/tracks/${Experimental_ID}"
    if [ -e "$outfile" ]
    then
        rm "$outfile"
    fi

    echo "Getting $outfile"
    ext=$(echo "$outfile" | cut -f2- -d'.')
    if [ "$ext" = "bw" ]
    then
        echo bigwig
        # BigWig Download URL:
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bw/${Experimental_ID} -O "$outfile"
    fi

    if [ "$ext" = "05.bb" ]
    then
        echo "bigbed"
        # Peak-call (BigBed) Download URL: #(Threshold = 05, 10, or 20)
        wget https://chip-atlas.dbcls.jp/data/${Genome}/eachData/bb05/${Experimental_ID} -O "$outfile"
    fi
done
