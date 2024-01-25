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
    parser.add_argument('--infile_tf_pos_prom'
        ,type=str
        ,default='results/mm10/tensor_TFsvd1_tf_pos_prom.npy'
        ,help="Prom x TF matrix of nr. of peaks per prom and TF")
    parser.add_argument('--outfile'
        ,type=str
        ,default='results/mm10/tensor_TFsvd1_avg_over_prom.npy'
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
    
    # get 1st comp. of chip signal for each TF
    if os.path.isfile(args.outfile)==False:
        # initilize tf x pos x prom tensor
        X = np.load(args.infile_tf_pos_prom)

        # average over promoters
        Peaks = np.zeros([N_tf,N_pos])
        Peaks_empty = np.zeros([N_tf,N_pos])
        for t,tf in enumerate(TFs):
            print(np.round(t/N_tf,3))
            idx_peaks = N_peaks.loc[:,tf]>0
            Peaks[t,:] = X[t,:,idx_peaks].mean(axis=0)
            Peaks_empty[t,:] = X[t,:,~idx_peaks].mean(axis=0)

        # Normalize
        Peaks = Peaks/Peaks.max(axis=1)[:,None]
        Peaks_empty = Peaks_empty/Peaks_empty.max(axis=1)[:,None]

        np.save(args.outfile,Peaks)
        np.save(args.outfile.replace('.npy','_empty.npy'),Peaks_empty)
    else:
        Peaks = np.load(args.outfile)
        Peaks_empty = np.load(args.outfile.replace('.npy','_empty.npy'))

    # remove high signal peaks from Peaks


    # plot mean chip signal per promoter for each tf
    if args.genome=='mm10':
        rows = 25
        cols = 17
    else:
        n = np.ceil(np.sqrt(N_tf)).astype(int)
        rows = n
        cols = n

    fig, axes = plt.subplots(rows,cols,sharex=True,sharey=True)

    x = np.linspace(-int(args.window_kb)*1000,int(args.window_kb)*1000,N_pos+1)
    x = .5*(x[:-1]+x[1:])
    for i,tf in enumerate(TFs):
        
        ax = axes.flatten()[i]

        ax.plot(x,Peaks[i,:],color='k',linewidth=0.5)
        # add point at tss
        ax.plot(0,np.mean(Peaks[i,[49,50]]),'r.',markersize=3)
        ax.set_title(f'{tf} ({np.sum(N_peaks.loc[:,tf]>0)})')
        # remove axis
        ax.axis('off')
    fig.set_size_inches([cols*4,rows*3])
    plt.tight_layout()
    fig.savefig(f'results/fig/{args.genome}/mean_chip_signal_per_promoter_with_peaks_each_tf.pdf')


    # plot mean chip signal per promoter for each tf
    fig, axes = plt.subplots(rows,cols,sharex=True,sharey=True)

    for i,tf in enumerate(TFs):
        
        ax = axes.flatten()[i]

        ax.plot(x,Peaks_empty[i,:],color='k',linewidth=0.5)
        ax.set_title(f'{tf} ({np.sum(N_peaks.loc[:,tf]==0)})')
        # remove axis
        ax.axis('off')
    fig.set_size_inches([cols*4,rows*3])
    plt.tight_layout()
    fig.savefig(f'results/fig/{args.genome}/mean_chip_signal_per_promoter_each_tf_no_peaks.pdf')