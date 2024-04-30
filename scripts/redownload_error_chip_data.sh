#!/bin/bash

genome=$1
infile=$2

mapfile -t to_download  < <(cat "$infile")

# pass if file is empty
if [ ${#to_download[@]} -eq 0 ]
then
    echo "No files to download"
    exit 0
else
    echo "Downloading ${#to_download[@]} files"

    # get files
    for Experimental_ID in "${to_download[@]}"
    do
        
        outfile="resources/tracks/${genome}/${Experimental_ID}"
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
            wget https://chip-atlas.dbcls.jp/data/"${genome}"/eachData/bw/"${Experimental_ID}" -O "$outfile"
        fi

        if [ "$ext" = "05.bb" ]
        then
            echo "bigbed"
            # Peak-call (BigBed) Download URL: #(Threshold = 05, 10, or 20)
            wget https://chip-atlas.dbcls.jp/data/"${genome}"/eachData/bb05/"${Experimental_ID}" -O "$outfile"
        fi
    done
    
fi