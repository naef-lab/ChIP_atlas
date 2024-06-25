import pandas as pd
import numpy as np
import argparse
import multiprocessing as mp
from functools import partial

def fill_in_peaks(p,promoterome,peak_table,N_pos):
    print(p)
    
    chr = promoterome.loc[p,'chr']
    start = promoterome.loc[p,'start']
    end = promoterome.loc[p,'end']
    peaks = peak_table[(peak_table['chr']==chr) & (peak_table['start']>=start) & (peak_table['end']<=end)]
    Peaks_p = np.zeros(N_pos)
    for i in peaks.index:
        start_pos = max(peaks.loc[i,'start'] - start,0)
        end_pos = min(peaks.loc[i,'end'] - start,N_pos)
        Peaks_p[start_pos:end_pos] += peaks.loc[i,'score']

    return Peaks_p
    

if __name__ == '__main__':
    
    genome = 'mm10'
    window_kb = 5
    infile_promoterome = f"/home/jbreda/Promoterome/results/{genome}/promoterome_pm{window_kb}kb_filtered.bed"
    promoterome = pd.read_csv(infile_promoterome,sep='\t')

    promoterome = promoterome.iloc[:24,:]

    infile_peak_table = f"/home/jbreda/ChIP_atlas/results/{genome}/peak_table_pm{window_kb}kb.txt"
    peak_table = pd.read_csv(infile_peak_table,sep='\t')

    # split peaks per prom
    Peaks_dic = {}
    for p in promoterome.index:
        print(p)
        chr = promoterome.loc[p,'chr']
        start = promoterome.loc[p,'start']
        end = promoterome.loc[p,'end']
        peaks = peak_table[(peak_table['chr']==chr) & (peak_table['start']>=start) & (peak_table['end']<=end)]
        peaks.start -= start
        peaks.end -= start
        Peaks_dic[p] = peaks

    N_prom = promoterome.shape[0]
    N_pos = 2*window_kb*1000

    # run loop in parralel
    threads = 24
    pool = mp.Pool(threads)

    outs = pool.map(partial(fill_in_peaks, promoterome=promoterome,peak_table=peak_table,N_pos=N_pos), promoterome.index)
    X = np.array(outs)

                  
    #X = np.array(out).reshape([N_tf,N_pos,N_prom])
    #for t in range(N_tf):
    #    # flip promoters on minus strand
    #    X[t,:,idx_minus_strand] = X[t,:,idx_minus_strand][:,::-1]
#

    #Peaks = np.zeros((N_prom,N_pos))

    #for p in promoterome.index:
    #    print(p)
#
    #    chr = promoterome.loc[p,'chr']
    #    start = promoterome.loc[p,'start']
    #    end = promoterome.loc[p,'end']
    #    peaks = peak_table[(peak_table['chr']==chr) & (peak_table['start']>=start) & (peak_table['end']<=end)]
#
    #    for i in peaks.index:
    #        start_pos = peaks.loc[i,'start'] - start
    #        end_pos = peaks.loc[i,'end'] - start
    #        Peaks[p,start_pos:end_pos] += peaks.loc[i,'score']
    #    
