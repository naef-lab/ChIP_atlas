import pandas as pd

configfile: 'config/chip_seq.yaml'


wildcard_constraints:
    genome="|".join(config['Genome']),
    version="|".join(config['Version'])


rule all:
    input:
        #expand("genome/{genome}/GeneName_Synonyms_dict.{ext}",genome=config['Genome'],ext=['txt','pkl']),
        expand("resources/{genome}/TF_list.txt",genome=config['Genome']),
        expand("resources/{genome}/TF_list_no_go_terms.txt",genome=config['Genome']),
        #expand("results/fig/hist_antigen_class_per_genome_{version}.pdf",version=['v2','v3']),
        #expand("results/fig/hist_experiment_QC_{version}.pdf",version=['v2','v3']),
        #expand("results/fig/hist_peaks_per_unique_mapped_read_{version}.pdf",version=['v2','v3']),
        #expand("resources/experimentList_{version}_{genome}_TFs_only_QC_filtered.tab",genome=['mm10','hg38'],version=['v2']),
        #expand("log/download_chip_data_{version}_{genome}.log",genome=['mm10','hg38'],version=['v3']),
        #expand("resources/to_redowload_{version}_{genome}.txt",genome=['mm10','hg38'],version=['v3']),
        #expand("log/redownload_error_chip_data_{version}_{genome}.log",genome=['mm10','hg38'],version=['v3'])

# Make GeneID to GeneName dictionary
rule GeneID_GeneName_Synonym_dict:
    input:
        "genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz"
    output:
        table="genome/{genome}/GeneName_Synonyms_dict.txt",
        dict="genome/{genome}/GeneName_Synonyms_dict.pkl"
    shell:
        "python scripts/get_GeneID_GeneName_Synonym_dict.py --infile {input} --outfile_table {output.table} --outfile_dict {output.dict}"

rule get_TF_list:
    input:
        geneid_genename_synonym="genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz",
        synonym_genename="genome/{genome}/GeneName_Synonyms_dict.pkl",
        TFs = lambda wildcards: config['Curated_TF_list'][wildcards.genome],
        GO_terms=lambda wildcards: expand("resources/GO_terms/{genome}/{go_term}.txt",go_term=config['GO_list'],genome=wildcards.genome)
    output:
        TFs="resources/{genome}/TF_list.txt"
    shell:
        """
        python scripts/get_TF_list.py --geneid_genename_synonym_table {input.geneid_genename_synonym} \
                                      --synonym_genename_dict {input.synonym_genename} \
                                      --infile_tf_list {input.TFs} \
                                      --infiles_GO_terms {input.GO_terms} \
                                      --genome {wildcards.genome} \
                                      --outfile {output.TFs}
        """

rule get_TF_list_no_go_terms:
    input:
        geneid_genename_synonym="genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz",
        synonym_genename="genome/{genome}/GeneName_Synonyms_dict.pkl",
        TFs = lambda wildcards: config['Curated_TF_list'][wildcards.genome],
    output:
        TFs="resources/{genome}/TF_list_no_go_terms.txt"
    shell:
        """
        python scripts/get_TF_list.py --geneid_genename_synonym_table {input.geneid_genename_synonym} \
                                      --synonym_genename_dict {input.synonym_genename} \
                                      --infile_tf_list {input.TFs} \
                                      --genome {wildcards.genome} \
                                      --outfile {output.TFs}
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
        tfs="resources/{genome}/TF_list_{version}.txt",
        geneid_genename_synonym="genome/{genome}/{genome}_ENSID_Genename_synonyms.txt.gz",
        synonym_genename="genome/{genome}/GeneName_Synonyms_dict.pkl",
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
        python scripts/get_experimentList_TF.py --infile_chip {input.chip} \
                                                --infile_tfs {input.tfs} \
                                                --geneid_genename_synonym_table {input.geneid_genename_synonym} \
                                                --synonym_genename_dict {input.synonym_genename} \
                                                --outfile {output.chip} \
                                                --genome {wildcards.genome} \
                                                --th_reads {params.th_reads} \
                                                --th_mapped_reads {params.th_mapped_reads} \
                                                --th_duplicates {params.th_duplicates} \
                                                --th_peaks {params.th_peaks} \
                                                --th_exp_per_tf {params.th_exp_per_tf}
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
        "resources/to_redowload_{version}_{genome}.txt"
    shell:
        """
        python scripts/check_tracks.py --genome {wildcards.genome} --infile {input} --outfile {output}
        """

# Redownload error tracks
rule redownload_error_chip_data:
    input:
        "resources/to_redowload_{version}_{genome}.txt"
    output:
        "log/redownload_error_chip_data_{version}_{genome}.log"
    shell:
        """
        bash scripts/redownload_error_chip_data.sh {wildcards.genome} {input} > {output}
        """