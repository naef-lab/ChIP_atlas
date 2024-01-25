import numpy as np
import pandas as pd
import os
import h5py
import matplotlib.pyplot as plt
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Make tensor of SVD 1st comp. of Chip signal')
    parser.add_argument('--infile_chip'
        ,type=str
        ,default="resources/experimentList_mm10_TFs_only_QC_filtered.tab"
        ,help="Chip experiment table")
    parser.add_argument('--infile_promoterome'
        ,type=str
        ,default="/home/jbreda/Promoterome/results/mm10/promoterome_pm2kb_filtered.bed"
        ,help="Promoterome bed file")
    parser.add_argument('--infile_tf_peak_per_prom'
        ,type=str
        ,default='results/mm10/nr_of_peaks_per_prom_per_tf.tsv'
        ,help="Prom x TF matrix of nr. of peaks per prom and TF")
    parser.add_argument('--infile_tensor'
        ,type=str
        ,default='results/mm10/tensor_TFsvd1_tf_pos_prom.npy'
        ,help="hdf5 with SVD 1st comp. tensor")
    parser.add_argument('--outfig'
        ,type=str
        ,default='results/fig/mm10/mean_chip_signal_per_promoter_each_tf.pdf'
        ,help="Figure with Chip signal distribution rel. to tss")
    parser.add_argument('--genome'
        ,type=str
        ,choices=['hg38','mm10']
        ,default='mm10'
        ,help="Genome version")
    parser.add_argument('--window_kb'
        ,type=int
        ,default=2
        ,help="Window size in kb")
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

    # get TF peaks per prom
    N_peaks = pd.read_csv(args.infile_tf_peak_per_prom,
                          sep='\t',
                          index_col=0)
    
    # get only promoters in promoterome
    prom_id = promoterome.gene + '_' + promoterome.id
    N_peaks = N_peaks.loc[prom_id,:]

    # get dims
    N_prom = promoterome.shape[0]
    N_pos = 100
    N_tf = len(TFs)

    # get chip signal
    X = np.load(args.infile_tensor)

    #average over TFs
    X = X.mean(axis=0)
    X_capped = np.clip(X,0,2)

    # plot chip signal distribution
    fig,axes = plt.subplots(1,2,figsize=(16,5))

    ax = axes[0]
    ax.hist(X.flatten(),bins=200)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_title('Chip signal distribution averaged over TFs')
    ax.set_xlabel('Chip signal')
    ax.set_ylabel('Frequency')


    ax = axes[1]
    ax.hist(X_capped.flatten(),bins=200)
    ax.set_yscale('log')
    ax.set_title('Capped at 10')
    ax.set_xlabel('Chip signal')
    ax.set_ylabel('Frequency')

    fig.savefig(f'results/fig/{args.genome}/chip_signal_distribution_tf_average.pdf')
    
