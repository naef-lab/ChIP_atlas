import numpy as np
import pandas as pd
import pyBigWig as bw
from multiprocessing import Pool
from functools import partial
import argparse


def argparser():
    parser = argparse.ArgumentParser(description='Get regulatory regions')
    parser.add_argument('--genome', type=str, default='mm10', help='Genome version')
    parser.add_argument('--regulatory_regions', type=str, default="results/mm10/regulatory_regions_pm5kb.bed", help='Regulatory regions bed file')
    parser.add_argument('--chip_experiments', type=str, default="resources/experimentList_v3_mm10_TFs_only_QC_filtered.tab", help='ChIP-seq experiment table')
    parser.add_argument('--out_matrix', type=str, default='results/mm10/N_sf_pm5kb.npy', help='Output file with N_sf matrix')
    parser.add_argument('--threads', type=int, default=24, help='Number of threads')

    return parser

def get_N_s(TFs,chip_experiment,genome,coord):
    chr,start,end = coord
    N_s = np.zeros(len(TFs))
    for f,tf in enumerate(TFs):
        n = 0
        for exp in chip_experiment[chip_experiment == tf].index:
            bbfile = bw.open(f'resources/tracks/{genome}/{exp}.05.bb')
            if chr not in bbfile.chroms():
                continue
            val = bbfile.entries(chr, start, end)
            if val != None:
                N_s[f] += np.sum([int(v[2]) for v in val])
                n += 1
            bbfile.close()
        if n > 0:
            N_s[f] /= n
        
    return N_s


if __name__ == '__main__':
    # get arguments
    parser = argparser()
    args = parser.parse_args()

    # load regulatory regions and chip experiments, get TFs
    Regulatory_regions_bed = pd.read_csv(args.regulatory_regions, sep='\t')
    chip_experiment = pd.read_csv(args.chip_experiments , sep='\t',index_col=0)
    chip_experiment = chip_experiment.loc[:,'antigen']
    TFs = chip_experiment.unique()

    # get reg. regions coordinates
    COORD = []
    for s in Regulatory_regions_bed.index:
        chr = Regulatory_regions_bed.loc[s,'chr']
        start = Regulatory_regions_bed.loc[s,'start']
        end = Regulatory_regions_bed.loc[s,'end']
        COORD.append( (chr,start,end) )

    # get N_s in parallel for all reg. regions s
    with Pool(processes=args.threads) as pool:
        OUT = pool.map(partial(get_N_s,TFs,chip_experiment,args.genome), COORD)

    # concatenate and save N_sf matrix
    N_sf = np.array(OUT)
    np.save(args.out_matrix,N_sf)