import torch
import pandas as pd
import argparse
import pyBigWig
import numpy as np
import h5py

def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP signals for tf in tensor (N_promoters x N_positions X N_experiments).')
    parser.add_argument('--tf'
        ,required=True
        ,type=str
        ,help="TF to get experiments")
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['mm10','mm39','hg19','hg38']
        ,help="genome version")
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,help="window size (in kb)")
    parser.add_argument('--bin_size'
        ,default=20
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
    parser.add_argument('--outfile_failed'
        ,required=True
        ,type=str
        ,help="Output text file with tensor")
    

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()
    infile_promoter = f"~/Datastructure/results/{args.genome}/promoterome_{args.window_kb}kb.gff"

    # load promoters
    promoter = pd.read_csv(infile_promoter ,sep='\t')
    promoter.chr = promoter.chr.apply(lambda x: 'chr'+x)
                           
    # get tensor dimention and initialize
    N_prom = promoter.shape[0]
    win = int(args.window_kb*1000)
    N_pos = int(win/args.bin_size)
    N_exp = len(args.infiles_tf)

    X = torch.zeros([N_prom,N_pos,N_exp])
    X[:] = torch.nan
    failed_bw = []
    n=0
    for exp in args.infiles_tf:
        try:
            with pyBigWig.open(exp) as bw:
                for p in range(N_prom):
                    [chr,start,end] = promoter.loc[p,['chr','start','end']]
                    X[p,:,n] = torch.from_numpy( np.array( bw.stats(chr, start, end, type="mean", nBins=N_pos) ).astype(float) )
                n+=1
        except:
            failed_bw.append(exp)
            print("failed " + exp)
    
    # write tensor in out h5py file after removing failed expeiment
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset(args.tf,data=X[:,:,:n])

    # write failed bw exp id
    with open(args.outfile_failed,'w') as fout:
        for exp in failed_bw:
            fout.write(f'{exp}\n')
