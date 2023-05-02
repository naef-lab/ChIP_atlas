import numpy as np
import pandas as pd
import h5py
import pickle
import os
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

if __name__ == '__main__':

    # Get all TFs
    infile='resources/experimentList_mm10_TFs.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()

    N_prom = 30114
    N_pos = 100
    N_tf = TFs.shape[0]
    N_pc = 10

    for tf in TFs:
        
        infiles={'U':f'results/svd/{tf}_U.npy','S':f'results/svd/{tf}_S.npy','Vh':f'results/svd/{tf}_Vh.npy','rho':f'results/corr/{tf}_rho.npy'}
        outfig=f'results/fig/svd/{tf}.pdf'

        if np.all([os.path.exists(f) for f in infiles.values()]):
            print(tf)
            try:
                S = np.load(infiles['S'])
            except:
                print(F'{tf} not yet available')
                continue
            try:
                rho = np.load(infiles['rho'])
            except:
                print(F'{tf} not yet available')
                continue
            #U = np.load(infiles['U'])
            try:
                Vh = np.load(infiles['Vh'])
            except:
                print(F'{tf} not yet available')
                continue

            # sample hierarchical clustering based on corr matrix
            n = rho.shape[0]
            triu_idx = np.triu_indices(n,1)
            corr_dist = 1 - rho[triu_idx[0],triu_idx[1]]
            linkage = hierarchy.linkage(corr_dist,optimal_ordering=True)

            # Get pcs, loadings and explained variance
            #X_svd = U @ np.diag(S)
            Loadings = np.square(Vh)
            explained_var = np.square(S)
            explained_var /= sum(explained_var)

            # plot results (exp. corr, exp. loadings, Variance explained per comp)
            fig, axes = plt.subplots(1,4)

            ax = axes[0]
            R = hierarchy.dendrogram(linkage, ax=ax, above_threshold_color='y',orientation='left',labels=None)
            ax.set_yticklabels([])
            #fig.colorbar(pos, ax=ax, shrink=0.5,location='top')
            #ax.set_xlabel('Experiments')
            #ax.set_ylabel('Experiments')

            # reorder
            rho = rho[R['leaves'],:]
            rho = rho[:,R['leaves']]
            Loadings = Loadings[:,R['leaves']]


            ax = axes[1]
            pos = ax.imshow(rho, cmap='RdBu_r', interpolation=None,vmin=-1,vmax=1,origin='lower')
            fig.colorbar(pos, ax=ax, shrink=0.5,location='right')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xlabel('Experiments')
            ax.set_ylabel('Experiments')

            ax = axes[2]
            pos = ax.imshow(Vh.T, cmap='RdBu_r', interpolation=None,origin='lower')
            fig.colorbar(pos, ax=ax,shrink=0.5, location='right')
            ax.set_yticklabels([])
            ax.set_xlabel('Components')
            ax.set_ylabel('Experiments')

            ax = axes[3]
            ax.plot(np.arange(explained_var.shape[0]),explained_var,'.-')
            ax.set_yscale('log')
            ax.set_xlabel('Component')
            ax.set_ylabel('Frac of variance')
            
            fig.set_size_inches([4*12,12])
            plt.tight_layout()
            fig.savefig(outfig)
            plt.close(fig)


    

