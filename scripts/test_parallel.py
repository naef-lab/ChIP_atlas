import numpy as np
import pandas as pd
import argparse
import pyBigWig
import h5py
import multiprocessing
from functools import partial
import time


def fill_experiment(exp,promoter,N_pos):
    N_prom = 50
    X = np.zeros([N_prom,N_pos])
    X[:] = np.nan
    try:
        with pyBigWig.open(exp) as bw:
            for p in range(N_prom):
                [chr,start,end] = promoter.loc[p,['chr','start','end']]
                X[p,:] = np.array( bw.stats(chr, start, end, type="mean", nBins=N_pos) ).astype(float)
    except:
        X[p,:] = -1

    return X

def get_chip_table(genome):
    infile=f"resources/experimentList_{genome}_TFs_only_QC_filtered.tab"
    experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3])
    experiment_tf.columns = ['id','antigen']
    return experiment_tf

def get_experiments(genome,tf):
    experiment_tf = get_chip_table(genome)
    IDs = list(experiment_tf.loc[experiment_tf.antigen==tf,'id'])
    files = [f"resources/tracks/{genome}/{id}.bw" for id in IDs]
    return files

if __name__ == '__main__':

    genome = 'hg38'
    tf = 'STAT2'

    infile_promoter = f"/home/jbreda/Promoterome/results/hg38/promoterome_pm2kb_filtered.bed"

    # load promoters
    promoter = pd.read_csv(infile_promoter ,sep='\t')
    promoter.chr = promoter.chr.apply(lambda x: 'chr'+x)
                           
    # get tensor dimention and initialize
    N_prom = 50#promoter.shape[0]
    win = int(2*1000)
    N_pos = 100
    my_experiments = get_experiments(genome,tf)
    N_exp = len(my_experiments)
    
    start = time.time()

    pool = multiprocessing.Pool(20)
    outs = pool.map(partial(fill_experiment, promoter=promoter, N_pos=N_pos), my_experiments)
    X = np.array(outs)
    #to_keep = np.zeros(N_exp).astype('bool')
    #for n in range(N_exp):
    #    X[:,:,n] = outs[n]

    print(time.time() - start)    

    start = time.time()

    X_ = np.zeros([N_prom,N_pos,N_exp])
    X_[:] = np.nan
    to_keep_ = np.zeros(N_exp).astype('bool')
    for i,exp in enumerate(my_experiments):
        X_[:,:,i] = fill_experiment(exp,promoter,N_pos)

    print(time.time() - start)    


