#!/bin/bash

promoterome=resources/mm10_promoters_v2.gff.gz
for f in results/bedfiles/*.bed
do  
    outfile="${f//.bed/_prom.bed}"
    echo $outfile

    bedtools intersect -a $f -b $promoterome | bedtools sort > $outfile
done