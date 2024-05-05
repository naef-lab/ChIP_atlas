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
    parser.add_argument('--infile_tf_pos_prom'
        ,type=str
        ,default='results/mm10/svd/Window_pm5kb_bin_size_10/tensor_TFsvd1_pos_prom.hdf5'
        ,help="Prom x TF matrix of nr. of peaks per prom and TF")
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
        ,default=5
        ,help="Window size in kb")
    parser.add_argument('--bin_size'
        ,type=int
        ,default=10
        ,help="Bin size in bp")
    
    return parser.parse_args()

if __name__ == '__main__':

    # args
    args = parse_argument()

    # Get all TFs
    experiment_tf = pd.read_csv(args.infile_chip,sep='\t',usecols=[0,3])
    experiment_tf.columns = ['id','antigen']
    TFs = list(experiment_tf.antigen.unique())

    # plot mean chip signal per promoter for each tf
    if args.genome=='mm10':
        rows = 28
        cols = 16
    elif args.genome=='hg38':
        rows = 32
        cols = 33

    fig, axes = plt.subplots(rows,cols,sharex=True,sharey=True)

    with h5py.File(args.infile_tf_pos_prom,'r') as hf:
        [N_tf,N_pos,N_prom] = hf['tf_pos_prom'].shape

        x = np.linspace(-int(args.window_kb)*1000,int(args.window_kb)*1000,N_pos+1)
        x = .5*(x[:-1]+x[1:])
        idx_0 = np.where( np.abs(x) == np.min(np.abs(x)) )[0]

        for i,tf in enumerate(TFs):

            vals = np.nanmean(hf['tf_pos_prom'][i],axis=1)
            
            ax = axes.flatten()[i]

            ax.plot(x,vals,color='k',linewidth=0.5)
            # add point at tss
            ax.plot(0,np.mean(vals[idx_0]),'r.',markersize=3)
            ax.set_title(f'{tf}')
            # remove axis
            ax.axis('off')
        fig.set_size_inches([cols*4,rows*3])
        plt.tight_layout()
        fig.savefig(f'results/fig/{args.genome}/mean_chip_signal_per_promoter_with_peaks_each_tf.pdf')
