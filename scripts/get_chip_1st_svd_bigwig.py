import pandas as pd
import numpy as np
import h5py
import pyBigWig
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Compute svd of Chip signal accross experiments')
    parser.add_argument('--infile_svd'
        ,required=True
        ,type=str
        ,help="Chip svd")
    parser.add_argument('--promoterome'
        ,required=True
        ,type=str
        ,help="promoterome bed file")
    parser.add_argument('--chrom_size'
        ,required=True
        ,type=str
        ,help="chromosome size file")
    parser.add_argument('--bin_size'
        ,default=10
        ,type=int
        ,help="bin size")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Bigwig file with 1st svd")
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_argument()

    # get promoterome
    promoterome = pd.read_csv(args.promoterome,sep='\t')
    N_prom = promoterome.shape[0]

    # get 1st svd
    with h5py.File(args.infile_svd, 'r') as hf:
        u = hf['U'][:,0]
        s = hf['S'][0]
        N_pos = hf['U'].shape[0]//N_prom

    # reshape 1st svd as (Prom x Pos) matrix
    chip = np.reshape(u*s,[N_prom,N_pos]).T

    # get chrom size and initialize chr_values
    chr_size = []
    chr_values = {}
    CHR = promoterome.chr.unique()
    with open(args.chrom_size,'r') as f:
        for l in f.readlines():
            l = l.strip()
            chr = 'chr'+l.split('\t')[0]
            size = int(l.split('\t')[1])
            if chr in CHR:
                chr_size.append((chr,size))
                chr_values[chr] = np.zeros(int(np.ceil(size/args.bin_size)))

    # fill in chr_values with chip signal
    for p in range(N_prom):

        # skip if all nan
        if np.all(np.isnan( chip[:,p] )):
            continue

        # get absolute and relative idx
        chr = promoterome.at[p,'chr']
        start_idx = int(np.ceil(promoterome.at[p,'start']/args.bin_size))
        rel_idx =  np.where( ~np.isnan(chip[:,p]) )[0]
        idx = start_idx + rel_idx

        # add values
        chr_values[chr][idx] = chip[rel_idx,p]
        # if 2 promoters are overlaping the same bin, the chip signal is the same for both

    # write bigwig
    with pyBigWig.open(args.outfile, 'w') as bw:

        # add headers
        bw.addHeader(chr_size)

        # fill in the big wig chr by chr
        for chr in chr_values.keys():
            print(chr)
            
            # get indices with non zero values
            idx = np.where( chr_values[chr] != 0 )[0]

            # get coordinates and values
            chrs = idx.shape[0]*[chr]
            starts = idx*args.bin_size
            ends = starts + args.bin_size
            values = chr_values[chr][idx]

            # add entries
            bw.addEntries(chrs, starts, ends, values)