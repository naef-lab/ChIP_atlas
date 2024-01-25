import h5py
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':

    Genome = ['hg38','mm10']

    fig, axes = plt.subplots(1,2,figsize=(2*6,4))

    for genome in Genome:

        TFs = [tf.split('.')[0] for tf in os.listdir(f'results/{genome}/TF_tensors') if tf.endswith('.hdf5')]

        Sum_counts = []

        for tf in TFs:
            print(tf)
            infile = f'results/{genome}/TF_tensors/{tf}.hdf5'
                
            with h5py.File(infile,'r') as hf:
                Sum_counts.extend(np.nansum(hf['chip_prom_pos_exp'][:],axis=0).sum(axis=0))

        ax = axes[Genome.index(genome)]

        ax.hist(np.log10(Sum_counts),bins=100)
        ax.set_xlabel(r'$log_{10}$ Total counts')
        ax.set_ylabel('Number of experiments')

        ax.set_title(f'{genome} - median: {np.median(Sum_counts):.0f}')
    
    plt.tight_layout()
    fig.savefig(f'results/fig/total_count_per_exp.pdf')
    plt.close()


