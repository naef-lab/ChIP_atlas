import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import pyBigWig as bw
import torch
import argparse

def get_atac_data(chrom,start,end,run_ids):

    X = np.zeros([10000,run_ids.shape[0]])  

    for i,run in enumerate(run_ids):
        infile = f'/home/jbreda/ATACseq/results/mapping/{run}_coverage.bw'
        f = bw.open(infile)
        X[:,i] = f.values(chrom.replace('chr',''),start,end)

    return X

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def argparser():
    parser = argparse.ArgumentParser(description='Get regulatory regions')
    parser.add_argument('--genome', type=str, default='mm10', help='Genome version')
    parser.add_argument('--window_kb', type=int, default=5, help='Window size in kb')
    parser.add_argument('--promoterome', type=str, default='~/Promoterome/results/mm10/promoterome_pm5kb_filtered_clustered_sorted.bed', help='Promoterome file')
    parser.add_argument('--atac_samples', type=str, default="/home/jbreda/ATACseq/resources/SraRunTable.txt", help='ATAC-seq data')
    parser.add_argument('--chip_peak_sparse_tensor', type=str, default="results/mm10/Peak_tensors/Window_pm5kb/PeakCenter_sparse_tensor_prom_tf_position.pt" ,help='Sparse tensor with ChIP-seq peak center')
    parser.add_argument('--out_bedfile', type=str, default='results/mm10/default_output.bed' ,help='Output file with regulatory regions bed file')
    parser.add_argument('--out_matrix', type=str, default='results/mm10/default_output.npy', help='Output file with regulatory regions in boolean matrix')

    return parser


if __name__ == '__main__':
    # get arguments and parameters
    parser = argparser()
    args = parser.parse_args()
    w = 100
    prominence = {'chip':1,'atac':2}

    # load promoterome
    promoterome = pd.read_csv(args.promoterome, sep='\t')

    # load atac sample names
    if args.genome == 'mm10':
        atac_samples = pd.read_csv(args.atac_samples,sep=',')
        atac_run_ids = atac_samples['Run']
        Tracks = ['chip','atac']
    else:
        Tracks = ['chip']

    # load chip peak sparse tensor
    Chip_peaks_prom_tf_pos = torch.load(args.chip_peak_sparse_tensor)

    # create reg. region matrix and output bed file
    Regulatory_regions = np.zeros([promoterome.shape[0],args.window_kb*2000]).astype(bool)
    Regulatory_regions_bed = pd.DataFrame(np.zeros([0,6]))
    Regulatory_regions_bed.columns = ['chr','start','end','name','score','strand']

    for prom in promoterome.index:
        if prom % 100 == 0:
            print(np.round(prom/len(promoterome.index)*100,2),'%')

        # get promoter coordinates
        chrom = promoterome.at[prom,'chr']
        start = promoterome.at[prom,'start']
        end = promoterome.at[prom,'end']
        strand = promoterome.at[prom,'strand']

        # initialize signals
        signal = {}
        peaks = {}
        peak_properties = {}

        # get chip peaks signal
        Chip_peaks_tf_pos = Chip_peaks_prom_tf_pos[prom].to_dense()
        if strand == '-':
            Chip_peaks_tf_pos = Chip_peaks_tf_pos.flip(0)
        Chip_peaks_tf_pos = Chip_peaks_tf_pos.numpy().sum(0)/1000
        signal['chip'] = moving_average(Chip_peaks_tf_pos, w)*w

        # get atac signal
        if args.genome == 'mm10':
            X = get_atac_data(chrom,start,end,atac_run_ids.values)
            signal['atac'] = np.max(X,axis=1)

        # find peaks in chip and atac signals
        my_peaks = np.zeros(args.window_kb*2000).astype(bool)

        for track in Tracks:
            peaks[track], peak_properties[track] = find_peaks(signal[track],width=20,distance=20,prominence=prominence[track],height=0)
            for l,r in zip( np.floor(peak_properties[track]['left_ips']).astype(int), np.ceil(peak_properties[track]['right_ips']).astype(int) ):
                my_peaks[l:r] = True

        # close small gaps
        peak_idx = np.where(my_peaks)[0]
        to_add = np.array([])
        for i,p in enumerate(peak_idx[:-1]):
            d = peak_idx[i+1] - p
            if (d > 1) and (d < 50):
                to_add = np.concatenate([to_add,np.arange(p+1,peak_idx[i+1])])
        my_peaks[to_add.astype(int)] = True
        peak_idx = np.where(my_peaks)[0]

        Regulatory_regions[prom,:] = my_peaks

        # add peaks to bed file
        reg_start = None
        reg_end = None
        name = promoterome.at[prom,'id']
        s = 0
        for i in peak_idx:
            
            # if first region starts
            if reg_start is None:
                reg_start = i
                reg_end = i
                continue

            # if region continues expand
            if i == reg_end + 1:
                reg_end = i
            
            # if region ends add to bed and start new region
            else:
                s += 1
                Regulatory_regions_bed = pd.concat([Regulatory_regions_bed,pd.DataFrame({'chr':chrom,'start':start+reg_start,'end':start+reg_end,'name':name,'score':reg_end-reg_start,'strand':strand},index=[s])],axis=0)
                reg_start = i
                reg_end = i
        
        # add last region if any
        if reg_start is not None:
            s += 1
            Regulatory_regions_bed = pd.concat([Regulatory_regions_bed,pd.DataFrame({'chr':chrom,'start':start+reg_start,'end':start+reg_end,'name':name,'score':reg_end-reg_start,'strand':strand},index=[s])],axis=0)
            
    # fix types
    Regulatory_regions_bed.chr = Regulatory_regions_bed.chr.astype(str)
    Regulatory_regions_bed.start = Regulatory_regions_bed.start.astype(int)
    Regulatory_regions_bed.end = Regulatory_regions_bed.end.astype(int)
    Regulatory_regions_bed.name = Regulatory_regions_bed.name.astype(str)
    Regulatory_regions_bed.score = Regulatory_regions_bed.score.astype(int)
    Regulatory_regions_bed.strand = Regulatory_regions_bed.strand.astype(str)

    # save output
    Regulatory_regions_bed.to_csv(args.out_bedfile,sep='\t',index=False,header=True)
    np.save(args.out_matrix,Regulatory_regions)