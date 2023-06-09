import numpy as np
import pandas as pd
import argparse
import pyBigWig
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

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()
    infile_promoter = f"/home/jbreda/Promoterome/results/{args.genome}/promoterome_pm{args.window_kb}kb_filtered.bed"

    # load promoters
    promoter = pd.read_csv(infile_promoter ,sep='\t')
    promoter.chr = promoter.chr.apply(lambda x: 'chr'+x)
                           
    # get tensor dimention and initialize
    N_prom = promoter.shape[0]
    win = int(args.window_kb*1000)
    N_pos = int(win/args.bin_size)
    N_exp = len(args.infiles_tf)

    X = np.zeros([N_prom,N_pos,N_exp])
    X[:] = np.nan
    exp_id = []
    failed_bw = []
    n=0
    for exp in args.infiles_tf:
        print(f'{exp}... ',end='')
        try:
            with pyBigWig.open(exp) as bw:
                for p in range(N_prom):
                    [chr,start,end] = promoter.loc[p,['chr','start','end']]
                    X[p,:,n] = np.array( bw.stats(chr, start, end, type="mean", nBins=N_pos) ).astype(float)
                n+=1
            exp_id.append(exp.split('/')[3].split('.')[0])
            print('done')
        except:
            failed_bw.append(exp.split('/')[3].split('.')[0])
            print('failed')
    
    X = X[:,:,:n]
    print(X.shape)
    # write tensor in out h5py file after removing failed expeiment
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset('chip_prom_pos_exp',data=X)
        hf.create_dataset('experiment_id',data=exp_id)
        if len(failed_bw)>0:
            hf.create_dataset('failed_experiment',data=failed_bw)
