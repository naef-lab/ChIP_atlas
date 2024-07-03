import numpy as np
import pandas as pd
import pyBigWig
from multiprocessing import Pool
from functools import partial
import torch

import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP Peaks for tf in a sparse tensor (N_promoters x N_TF x N_positions).')
    parser.add_argument('--chip_experiment_table'
        ,required=True
        ,type=str
        ,help="ChIP-seq experiment table")
    parser.add_argument('--infile_promoterome'
        ,default="/home/jbreda/Promoterome/results/mm10/promoterome_pm1kb_filtered.bed"
        ,type=str
        ,help="Promoterome bed files")
    parser.add_argument('--outfile_sparse_tensor'
        ,required=True
        ,type=str
        ,help="Output sparse tensor file")
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,choices=[1,2,3,5]
        ,help="window size (in kb)")
    parser.add_argument('--threads'
        ,default=1
        ,type=int
        ,help="Number of threads")

    return parser.parse_args()


# get the peak positions of the TFs in the promoters
def find_tf_peaks_in_promoter(promoterome,experiment_tf,window_kb,tf_j):
    tf = tf_j[0]
    j = tf_j[1]

    # get TF experiments IDs
    IDs = list(experiment_tf.loc[experiment_tf.antigen==tf,'id'])

    # initialize IJKval matrix
    IJKval = np.zeros([4,0]) # n_dim (i,j,k,val) x n_peaks with i=promoter, j=tf, k=position, val=score

    # loop over all promoters
    for i,prom in enumerate( promoterome.index ):
        # get promoter coordinates
        [chr,start,end,strand] = promoterome.loc[prom,['chr','start','end','strand']].values

        # initialize score matrix (n_experiments x window size)
        score = np.zeros([len(IDs),window_kb*2*1000])*np.nan

        # loop over all experiments
        for l,id in enumerate(IDs):

            infile = f"resources/tracks/mm10/{id}.05.bb"
            # open the bigwig file (continue if it doesn't open or exists)
            try:
                bb = pyBigWig.open(infile)
            except:
                continue
        
            # check if the chromosome is in the bigwig file and get the entries
            if not chr in bb.chroms():
                continue

            # get the entries in the window size around the promoter exit if empty
            entries = bb.entries(chr, start, min(end,bb.chroms(chr)))
            if entries == None:
                continue

            # fill in the score matrix with the score of the experiment at the peak relative positions
            for e in entries:
                score[l,(e[0]-start):(e[1]-start)] = int(e[2])
            
            # close the bigwig file
            bb.close()

        # if there are no peaks for this TF in this promoter (or bb file couldn't open), then skip
        if np.isnan(score).all():
            continue

        # flip the score matrix if the promoter is on the negative strand
        if strand == '-':
            score = score[:,::-1]

        # find the peak positions in the score matrix
        peak_idx = np.where(np.any(~np.isnan(score),0))[0]
        # find the start and end of the peaks
        dp = np.diff(peak_idx)
        dp = np.insert(dp,0,0)
        dp = np.append(dp,0)
        ends = peak_idx[dp[1:] != 1]
        starts = peak_idx[dp[:-1] != 1]
        peaks = np.concatenate([starts.reshape(-1,1),ends.reshape(-1,1)],1)
        # add 1 to ends to make it half open
        peaks[:,1] += 1
        # get the score and position and save in Peaks matrix
        for peak in peaks:
            mean_score = np.nanmean(score[:,peak[0]:peak[1]])
            k = int((peak[0]+peak[1])/2)
            ijkval = np.array([i,j,k,mean_score])[:,None]
            IJKval = np.concatenate([IJKval,ijkval],1)

        
    return IJKval


if __name__ == '__main__':
    # parameters
    args = parse_argument()

    # load promoters
    promoterome = pd.read_csv(args.infile_promoterome ,sep='\t')
    CHR = promoterome.chr.unique()
    STRAND = ['+','-']

    # load chip experiments table
    experiment_tf = pd.read_csv(args.chip_experiment_table,sep='\t',usecols=[0,3])
    experiment_tf.columns = ['id','antigen']

    # get the TFs and their index
    TFs = experiment_tf.antigen.unique()
    TF_IDX = np.concatenate([TFs[:,None],np.arange(len(TFs))[:,None]],1)

    # get the peak positions of the TFs in the promoters run in parallel for all TFs
    with Pool(processes=args.threads) as pool:
        OUT = pool.map(partial(find_tf_peaks_in_promoter,promoterome,experiment_tf,args.window_kb), TF_IDX)
    
    # concatenate the results in one matrix (i,j,k,val) x N_peaks, where i=promoter, j=tf, k=position, val=score
    IJKval = np.zeros([4,0])
    for i,out in enumerate(OUT):
        IJKval = np.concatenate([IJKval,out],1)

    # make sparse tensor
    indices = IJKval[:3].astype(int)
    values = IJKval[3]
    my_size = [promoterome.shape[0],TFs.shape[0],args.window_kb*2*1000]
    PromTfPos = torch.sparse_coo_tensor(indices, values, size=my_size)

    # save tensor
    torch.save(PromTfPos, args.outfile_sparse_tensor)


