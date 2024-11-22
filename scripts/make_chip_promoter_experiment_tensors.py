import numpy as np
import pandas as pd
import argparse
import pyBigWig
import h5py
import multiprocessing
from functools import partial

def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP signals for tf in tensor (N_promoters x N_positions X N_experiments).')
    parser.add_argument('--promoterome'
        ,required=True,
        type=str,
        help="Promoterome file")
    parser.add_argument('--threads'
        ,required=True
        ,type=int
        ,help="nr. of threads to use")
    parser.add_argument('--tf'
        ,required=True
        ,type=str
        ,help="TF to get experiments")
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['mm10','hg38']
        ,help="genome version")
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,help="window size (in kb)")
    parser.add_argument('--bin_size'
        ,default=10
        ,type=int
        ,help="bin size")
    parser.add_argument('--infiles_tf'
        ,required=True
        ,type=str
        ,nargs='+'
        ,help="experiments coverage files")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output hfd5 file with tensor")

    return parser.parse_args()

def fill_experiment(exp,promoter,N_pos):
    N_prom = promoter.shape[0]
    X = np.zeros([N_prom,N_pos])
    try:
        with pyBigWig.open(exp) as bw:
            for p in range(N_prom):
                [chr,start,end] = promoter.loc[p,['chr','start','end']]
                if chr in bw.chroms():
                    X[p,:] = np.array( bw.stats(chr, start, end, type="mean", nBins=N_pos) ).astype(float)
                else:
                    X[p,:] = -1
    except:
        X[:] = -1

    return X

if __name__ == '__main__':

    args = parse_argument()

    # load promoters
    promoter = pd.read_csv(args.promoterome ,sep='\t')

    # get tensor dimention and initialize
    N_prom = promoter.shape[0]
    win_size = promoter.at[0,'end'] - promoter.at[0,'start']
    N_pos = int(win_size/args.bin_size)
    N_exp = len(args.infiles_tf)

    # fill tensor
    pool = multiprocessing.Pool(args.threads)
    outs = pool.map(partial(fill_experiment, promoter=promoter, N_pos=N_pos), args.infiles_tf)

    # remove failed experiments
    X = np.array(outs)
    X = np.transpose(X,axes=(1,2,0))
    to_keep = np.nanmean(np.nanmean(X,axis=0),axis=0) > 0
    X = X[:,:,to_keep]
    
    # get experiment IDs
    IDs = np.array( [exp.split('/')[3].split('.')[0] for exp in args.infiles_tf] )
    exp_id = [str(e) for e in IDs[to_keep]]
    failed_bw = [str(e) for e in IDs[~to_keep]]

    # write tensor in out h5py file after removing failed expeiment
    with h5py.File(args.outfile,'w') as hf:
        d = hf.create_dataset('chip_prom_pos_exp',data=X)
        d.attrs['exp_id'] = ','.join(IDs)
        d.attrs['failed_bw'] = ','.join(failed_bw)
        d.attrs['TF'] = args.tf
        d.attrs['genome'] = args.genome
        d.attrs['window_kb'] = args.window_kb
        d.attrs['bin_size'] = args.bin_size
        d.attrs['promoterome'] = args.promoterome
        
        
