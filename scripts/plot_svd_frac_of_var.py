import numpy as np
import pandas as pd
import h5py
import pickle
import os
import matplotlib.pyplot as plt
from scipy import signal

if __name__ == '__main__':

    # Get all TFs
    infile='resources/experimentList_mm10_TFs.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()

    N_prom = 30114
    N_pos = 100
    N_tf = TFs.shape[0]
    N_pc = 10

    Frac_var = {}

    for tf in TFs:
        
        infiles={'U':f'results/svd/{tf}_U.npy','S':f'results/svd/{tf}_S.npy','Vh':f'results/svd/{tf}_Vh.npy','rho':f'results/corr/{tf}_rho.npy'}
        outfig=f'results/fig/svd/{tf}.pdf'

        if os.path.exists(infiles['S']):
            print(tf)
            try:
                S = np.load(infiles['S'])
            except:
                print(F'{tf} not yet available')
                continue
            explained_var = np.square(S)
            Frac_var[tf] = explained_var/explained_var[0]


    fig = plt.figure()
    ax = fig.add_subplot(131)
    for tf in Frac_var:
        f = Frac_var[tf]
        if f.shape[0]>5:
            #f[f<1] = 1
            x = np.linspace(0,1,len(f))
            #x = np.arange(1,len(f)+1)
            ax.plot(x,f)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    #ax.set_xlim([1,100])
    ax.set_ylim([1e-8,1])
    ax.set_xlabel('Component')
    ax.set_ylabel('Frac of variance')

    ax = fig.add_subplot(132)
    for tf in Frac_var:
        f = Frac_var[tf]
        f = f[:-1]
        if f.shape[0]>5:
            lfc = np.log2(f[:-1]/f[1:])
            lfc /= lfc[0]
            #lfc[lfc>1] = 1
            lfc[lfc<0.1] = 0
            w = 3
            y = np.convolve(lfc, np.ones(w), 'full') / w
            x = np.arange(y.shape[0])
            x = np.linspace(0,1,y.shape[0])
            ax.plot(x,y)
    ax.set_xlabel('Component')
    ax.set_ylabel('log2 FC pc var')

    ax = fig.add_subplot(133)
    #for tf in Frac_var:
    #f = Frac_var[tf]
    #lfc = np.log(f[:-1]/f[1:])
    #dlfc = lfc[1:] - lfc[:-1]
    #ax.plot(np.arange(dlfc.shape[0])+1,dlfc)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlabel('Component')
    ax.set_ylabel('d lfc')

    fig.set_size_inches([3*12,8])
    plt.tight_layout()
    fig.savefig(f'results/fig/svd/Frac_var.pdf')
    plt.close(fig)



        

