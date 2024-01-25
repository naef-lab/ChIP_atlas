import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool

def parse_argument():
    parser = argparse.ArgumentParser(description='Make Prom x TF matrix of nr. of peaks per')
    parser.add_argument('--infile_chip'
        ,dest='infile_chip'
        ,type=str
        ,default="resources/experimentList_mm10_TFs_only_QC_filtered.tab"
        ,help="Chip experiment table")
    parser.add_argument('--infile_promoterome'
        ,type=str
        ,default="/home/jbreda/Promoterome/results/mm10/promoterome_pm2kb_filtered.bed"
        ,help="Promoterome bed file")
    parser.add_argument('--outfile'
        ,type=str
        ,default='results/mm10/nr_peaks_per_prom_per_tf.tsv'
        ,help="output file: Prom x TF matrix of nr. of peaks per prom and TF")
    
    return parser.parse_args()

if __name__ == '__main__':

    # args
    args = parse_argument()
    
    promoterome = pd.read_csv(args.infile_promoterome,sep='\t')
    idx_minus_strand = (promoterome.strand=='-').values

    # Get all TFs
    experiment_tf = pd.read_csv(args.infile_chip,sep='\t',usecols=[0,3])
    experiment_tf.columns = ['id','antigen']
    TFs = list(experiment_tf.antigen.unique())

    # get dims
    N_prom = promoterome.shape[0]
    N_tf = len(TFs)


    N = np.zeros([N_prom,N_tf])

    for tf in TF:
        print(tf)