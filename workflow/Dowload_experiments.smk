import pandas as pd

configfile: 'config/chip_seq.yaml'


wildcard_constraints:
    genome="|".join(config['Genome']),
    version="|".join(config['Version'])


rule all:
    input:
        #expand("results/fig/hist_antigen_class_per_genome_{version}.pdf",version=['v2','v3']),
        #expand("results/fig/hist_experiment_QC_{version}.pdf",version=['v2','v3']),
        #expand("results/fig/hist_peaks_per_unique_mapped_read_{version}.pdf",version=['v2','v3']),
        #expand("resources/experimentList_{version}_{genome}_TFs_only_QC_filtered.tab",genome=['mm10','hg38'],version=['v2','v3']),
        expand("log/download_chip_data_{version}_{genome}.log",genome=['mm10','hg38'],version=['v2','v3']),
        expand("resources/to_dowload_{version}_{genome}.txt",genome=['mm10','hg38'],version=['v2','v3']),
        expand("log/redownload_error_chip_data_{version}_{genome}.log",genome=['mm10','hg38'],version=['v2','v3'])

# Make GeneID to GeneName dictionary
rule GeneID_GeneName_Synonym_dict:
    input:
        "genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz"
    output:
        "genome/{genome}/{genome}_GeneID_GeneName_Synonyms_dict.txt"
    shell:
        "python scripts/get_GeneID_GeneName_Synonym_dict.py --infile {input} --outfile {output}"

rule get_mm10_TF_list:
    input:
        TFs="genome/mm10/mm10_TF_ID_list.csv",
        gene_dict="genome/mm10/mm10_GeneID_GeneName_Synonyms_dict.txt"
    output:
        TFs="genome/mm10_TF_list.csv"
    shell:
        """
        python scripts/get_mm_TF_list.py --infile_tf {input.TFs} --infile_gene_dict {input.gene_dict} --outfile {outpu.TFs}
        cp {output} /bigdata/jbreda/genome/mm39_TF_list.csv
        """

# Plot all QC stats for all genomes
rule plot_chip_experiments_stats_QC:
    input:
        chip="resources/experimentList_{version}.tab"
    output:
        antigen_class_per_genome="results/fig/hist_antigen_class_per_genome_{version}.pdf",
        QC="results/fig/hist_experiment_QC_{version}.pdf",
        peaks_per_unique_mapped_reads="results/fig/hist_peaks_per_unique_mapped_read_{version}.pdf",
    params:    
        genomes=config['Genome'],
        th_n_exp_per_tf=config['Threshold']['n_exp_per_tf'],
        th_reads=config['Threshold']['n_reads'],
        th_mapped_reads=config['Threshold']['f_mapped_reads'],
        th_duplicates=config['Threshold']['f_duplicates'],
        th_peaks=config['Threshold']['n_peaks'],
    shell:
        """
        python scripts/plot_chip_atlas_stats.py --infile {input.chip} --outfig_antigen_class_per_genome {output.antigen_class_per_genome} --outfig_QC {output.QC} --outfig_peaks_per_unique_mapped_reads {output.peaks_per_unique_mapped_reads}  --genomes {params.genomes}  --th_reads {params.th_reads} --th_mapped_reads {params.th_mapped_reads} --th_duplicates {params.th_duplicates} --th_peaks {params.th_peaks}
        """

# Fliter experiment list based on QC thresholds
rule get_ChIP_experiments_TF_QC_filtered:
    input:
        chip="resources/experimentList_{version}.tab",
        tfs="genome/{genome}/{genome}_TF_list.csv",
        gene_dict="genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz"
    output:
        chip="resources/experimentList_{version}_{genome}_TFs_only_QC_filtered.tab"
    params:
        th_reads=config['Threshold']['n_reads'],
        th_mapped_reads=config['Threshold']['f_mapped_reads'],
        th_duplicates=config['Threshold']['f_duplicates'],
        th_peaks=config['Threshold']['n_peaks'],
        th_exp_per_tf=config['Threshold']['n_exp_per_tf']
    shell:
        """
        python scripts/get_experimentList_TF.py --infile_chip {input.chip} --infile_tfs {input.tfs} --infile_gene_dict {input.gene_dict} --outfile {output.chip} --genome {wildcards.genome} --th_reads {params.th_reads} --th_mapped_reads {params.th_mapped_reads} --th_duplicates {params.th_duplicates} --th_peaks {params.th_peaks} --th_exp_per_tf {params.th_exp_per_tf}
        """

# First download data: download_chip_data.sh -> check_tracks.py -> redownload_error_chip_data.sh
rule download_chip_data:
    input:
        "resources/experimentList_{version}_{genome}_TFs_only_QC_filtered.tab"
    output:
        "log/download_chip_data_{version}_{genome}.log"
    shell:
        """
        bash scripts/download_chip_data.sh {wildcards.genome} {input} > {output}
        """

# Check if all tracks are downloaded
rule check_tracks:
    input:
        "resources/experimentList_{version}_{genome}_TFs_only_QC_filtered.tab"
    output:
        "resources/to_dowload_{version}_{genome}.txt"
    shell:
        """
        python scripts/check_tracks.py --genome {wildcards.genome} --infile {input} --outfile {output}
        """

# Redownload error tracks
rule redownload_error_chip_data:
    input:
        "resources/to_dowload_{version}_{genome}.txt"
    output:
        "log/redownload_error_chip_data_{version}_{genome}.log"
    shell:
        """
        bash scripts/redownload_error_chip_data.sh {wildcards.genome} {input} > {output}
        """