import pandas as pd
import numpy as np
import pyBigWig as bw
import argparse


def parse_argument():
    parser = argparse.ArgumentParser(description='get chromosome size file from bigwig')
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['mm10','hg38']
        ,help="genome")
    parser.add_argument('--chip_experiment'
        ,required=True
        ,type=str
        ,help="chip experiment file")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="chromosoze size file")
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_argument()

    chip = pd.read_csv(args.chip_experiment, sep='\t',index_col=0)

    # pick 1000 random experiments and get chrom size
    chr_size = pd.DataFrame()
    for id in chip.sample(1000).index:
        with bw.open(f'resources/tracks/{args.genome}/{id}.bw') as fin:
            chr_size = pd.concat([chr_size,pd.DataFrame.from_dict(fin.chroms(),orient='index',columns=[id])],axis=1)

    # check if all experiments have the same chrom size and write to file
    with open(args.outfile,'w') as fout:
        for chr in chr_size.index:
            vals = chr_size.loc[chr].values
            val = np.unique( vals[~np.isnan(vals)].astype(int) )

            if len(val) > 1:
                print(f'{chr} has multiple sizes: {val}, take max')
                val = val.max()
            
            fout.write(f'{chr}\t{val[0]}\n')
        

