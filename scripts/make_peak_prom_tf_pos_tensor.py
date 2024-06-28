import numpy as np
import pandas as pd
import pyBigWig
from multiprocessing import Pool
from functools import partial
import torch

import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP Peaks for tf in a sparse tensor (N_promoters x N_TF x N_positions).')
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['mm10','hg38']
        ,help="genome version")
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,choices=[1,2,3,5]
        ,help="window size (in kb)")
    parser.add_argument('--infile_promoterome'
        ,default="/home/jbreda/Promoterome/results/mm10/promoterome_pm1kb_filtered.bed"
        ,type=str
        ,help="Promoterome bed files")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output sparse tensor")

    return parser.parse_args()

if __name__ == '__main__':
    # parameters
    args = parse_argument()

    # load promoters
    promoterome = pd.read_csv(args.infile_promoter ,sep='\t')
    CHR = promoterome.chr.unique()

    # load chip experiments table
    infile=f"../resources/experimentList_{args.genome}_TFs_only_QC_filtered.tab"
    experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3])
    experiment_tf.columns = ['id','antigen']

    TFs = experiment_tf.antigen.unique()