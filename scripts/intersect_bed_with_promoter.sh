#!/bin/bash

promoterome=resources/mm10_promoters_v2_pm1kb.gff
for f in results/TF_bedfiles/*.bed
do  
    outfile="${f//TF_bedfiles/TF_prom_bedfiles}"
    echo $outfile

    bedtools intersect -a $f -b $promoterome | bedtools sort | uniq > $outfile
done