#!/bin/bash

# Get TF list and count
mapfile -t TF < <(find results/bedfiles/*.bed | cut -f1 -d'.' | cut -f3 -d'/')
N=${#TF[@]}

# prepare output table
outfile='results/Jaccard_peak_pairs.txt'
if [ -e $outfile ]
then
    rm $outfile
fi
touch $outfile

# write header 
for j in $(seq "$N")
do
    printf "\t%s" "${TF[$j-1]}" >> $outfile
done

# fill table with jaccard indices
for i in $(seq "$N")
do
    echo "$i"/"$N"

    # Write row name and fill row
    printf "\n%s" "${TF[$i-1]}" >> $outfile
    for j in $(seq "$N")
    do
        if [[ $j > $i ]]
        then 
            printf "\t%f" "$(bedtools jaccard -a results/bedfiles/"${TF[$i-1]}".bed -b results/bedfiles/"${TF[$j-1]}".bed | tail -1 | cut -f3)"  >> $outfile
        else
            printf "\t0" >> $outfile
        fi
    done
done

# last EOL
printf "\n" >> $outfile