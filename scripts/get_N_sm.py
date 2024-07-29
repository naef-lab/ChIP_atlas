import numpy as np
import pandas as pd
import h5py
import argparse

def argparser():
    parser = argparse.ArgumentParser(description='Get regulatory regions')
    parser.add_argument('--regulatory_regions', type=str, default="results/mm10/regulatory_regions_pm5kb.bed", help='Regulatory regions bed file')
    parser.add_argument('--promoterome', type=str, default='~/Promoterome/results/mm10/promoterome_pm5kb_filtered_clustered_sorted.bed', help='Promoterome bed file')
    parser.add_argument('--jaspar_conv_tensor', type=str, default="~/Jaspar/results/mm10/Window_pm5kb/convolution_PromSeq_PWM.hdf5", help='Convolution tensor with Jaspar PWMs in h5')
    parser.add_argument('--threshold', type=float, default=0.5, help='Threshold for PWMs score')
    parser.add_argument('--out_matrix', type=str, default='results/mm10/N_sm_pm5kb.npy', help='Output file with N_sm matrix')

    return parser


if __name__ == '__main__':
    # get arguments
    parser = argparser()
    args = parser.parse_args()

    # load regulatory regions and chip experiments, get TFs
    Regulatory_regions_bed = pd.read_csv(args.regulatory_regions, sep='\t')

    # load Promoterome
    Promoterome = pd.read_csv(args.promoterome, sep='\t')

    # add promoter index and gene name to Regulatory regions
    prom2idx =  dict(zip(Promoterome.id.values, Promoterome.index.values))
    prom2gene = dict(zip(Promoterome.id.values, Promoterome.gene.values))
    prom2prom_start = dict(zip(Promoterome.id.values, Promoterome.start.values))
    prom2prom_end = dict(zip(Promoterome.id.values, Promoterome.end.values))
    Regulatory_regions_bed['prom_idx'] = Regulatory_regions_bed['name'].map(prom2idx)
    Regulatory_regions_bed['gene'] = Regulatory_regions_bed['name'].map(prom2gene)
    Regulatory_regions_bed['prom_start'] = Regulatory_regions_bed['name'].map(prom2prom_start)
    Regulatory_regions_bed['prom_end'] = Regulatory_regions_bed['name'].map(prom2prom_end)

    # load motif matrix
    Jaspar = h5py.File(args.jaspar_conv_tensor, 'r')

    # get dimentions
    N_reg_regions = Regulatory_regions_bed.shape[0]
    N_motifs = Jaspar['convolution'].shape[1]
    N_sm = np.zeros((N_reg_regions,N_motifs))
    for s in Regulatory_regions_bed.index:
        if s%1000 == 0:
            print(s/N_reg_regions*100, '%')

        # get coordinates
        p = Regulatory_regions_bed.loc[s,'prom_idx']
        start = Regulatory_regions_bed.loc[s,'start'] - Regulatory_regions_bed.loc[s,'prom_start']
        end = Regulatory_regions_bed.loc[s,'end'] - Regulatory_regions_bed.loc[s,'prom_start']
        
        # get motif matrix, apply threshold and sum
        X = Jaspar['convolution'][p,:,start:end]
        X[X<args.threshold] = 0
        N_sm[s] = X.sum(axis=1)

    np.save(args.out_matrix,N_sm)