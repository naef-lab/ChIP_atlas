#!/bin/bash

# params:
tf="Ar"
genome="mm10"
window_kb=5
bin_size=10

# input:
promoterome="/home/jbreda/Promoterome/results/${genome}/promoterome_pm${window_kb}kb_filtered_clustered_sorted.bed"

exp_list="resources/experimentList_Ar.tab"
# make a list of input files names "resources/tracks/${genome}/${id}.bw" with id in IDs
infiles_tf=$(cut -f 1 ${exp_list} | tr ' ' '\n' | sed "s|^|resources/tracks/${genome}/|g" | sed "s|$|.bw|g" | tr '\n' ' ')

#exp_list="resources/experimentList_v3_${genome}_TFs_only_QC_filtered.tab"

# output:
tensor="results/${genome}/Chip_tensors/Window_pm${window_kb}kb_bin_size_${bin_size}/${tf}.hdf5"


# run the script
threads=12
python scripts/make_chip_promoter_experiment_tensors.py --promoterome ${promoterome} \
                                                        --threads ${threads} \
                                                        --tf ${tf} \
                                                        --genome ${genome} \
                                                        --window_kb ${window_kb} \
                                                        --bin_size ${bin_size} \
                                                        --outfile ${tensor} \
                                                        --infiles_tf "${infiles_tf}"



                                                        