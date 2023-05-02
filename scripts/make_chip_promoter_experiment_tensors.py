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
    parser.add_argument('--infile_promoters'
        ,required=True
        ,type=str
        ,help="file with intronic reads per barcode")
    parser.add_argument('--bin_size'
        ,default=20
        ,type=int
        ,help="bin size")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output hfd5 file with tensor")
    parser.add_argument('--infiles_tf'
        ,required=True
        ,type=str
        ,nargs='+'
        ,help="experiments coverage files")

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()

    # load promoters
    promoter = pd.read_csv(args.infile_promoters ,sep=r'\t|\;',header=None,usecols=[0,3,4,6,8,9,10,11,12],engine='python')
    promoter.columns = ['chr','start','end','strand','ID','Name','Members','Annotations','CpG_class']
    for col in promoter.columns[4:]:
        promoter.loc[:,col].replace(to_replace='\"', value='', regex=True, inplace=True)
        promoter.loc[:,col].replace(to_replace=f'{col}\=',value='',regex=True, inplace=True)

    # get tensor dimention and initialize
    N_prom = promoter.shape[0]
    N_pos = int(2000/args.bin_size)
    N_exp = len(args.infiles_tf)

    X = torch.zeros([N_prom,N_pos,N_exp])
    X[:] = torch.nan
    failed_bw = []
    for k,exp in enumerate(args.infiles_tf):
        try:
            with pyBigWig.open(exp) as bw:
                for p in range(N_prom):
                    [chr,start,end] = promoter.loc[p,['chr','start','end']]
                    X[p,:,k] = torch.from_numpy( np.array( bw.stats(chr, start, end, type="mean", nBins=N_pos) ).astype(float) )
        except:
            print(f"Unable to open {exp}")
            failed_bw.append(k)

    # write tensor in out h5py file
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset(args.tf,data=X)

    # write failed bw rows and files
    with open(f'results/TF_tensors/{args.tf}_failed.txt','w') as fout:
        for k in failed_bw:
            fout.write(f'{k}\t{args.infiles_tf[k]}\n')
