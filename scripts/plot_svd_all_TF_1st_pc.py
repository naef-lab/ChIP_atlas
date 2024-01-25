import numpy as np
import pandas as pd
import h5py
import os
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

if __name__ == '__main__':

    outfig = 'results/fig/svd_all_TFs_1st_comp.pdf'

    # Get all TFs
    infile='resources/experimentList_mm10_TFs.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()
    TFs[TFs!='Ctcf']

    N_prom = 30114
    N_pos = 100
    N_tf = TFs.shape[0]
    N_pc = 3

    X_tf = np.zeros([N_prom*N_pos,N_tf])
    

    for t,tf in enumerate(TFs):
        print(np.round(t/N_tf,3))
        
        infiles_svd={'U':f'results/svd/{tf}_U.npy','S':f'results/svd/{tf}_S.npy','Vh':f'results/svd/{tf}_Vh.npy','rho':f'results/corr/{tf}_rho.npy'}
        infile_X=f'results/TF_tensors/{tf}.hdf5'

        if np.all([os.path.exists(f) for f in infiles_svd.values()]):
            U = np.load(infiles_svd['U'])
            S = np.load(infiles_svd['S'])
            X_tf[:,t] = U[:,0] * S[0]

        else:
            with h5py.File(infile_X,'r') as hf:
                X_tf[:,t] = hf[tf][:][:,:,0].reshape([N_prom*N_pos])

    # remove TFs with nans everywhere
    idx_out = np.isnan(X_tf).all(axis=0)
    X_tf = X_tf[:,~idx_out]

    # replace nans with 0s and compute correlation matrix
    X_tf[np.isnan(X_tf)] = 0
    rho = np.corrcoef(X_tf.T)

    # TF hierarchical clustering based on corr matrix
    n = rho.shape[0]
    triu_idx = np.triu_indices(n,1)
    corr_dist = 1 - rho[triu_idx[0],triu_idx[1]]
    linkage = hierarchy.linkage(corr_dist,optimal_ordering=True)

    # plot results (exp. corr, exp. loadings, Variance explained per comp)
    fig, axes = plt.subplots(1,2)

    ax = axes[0]
    R = hierarchy.dendrogram(linkage, ax=ax, above_threshold_color='y',orientation='left',labels=None)
    ax.set_yticklabels([])

    # reorder corr matrix
    rho = rho[R['leaves'],:]
    rho = rho[:,R['leaves']]

    ax = axes[1]
    pos = ax.imshow(rho, cmap='RdBu_r', interpolation=None,vmin=-1,vmax=1,origin='lower')
    fig.colorbar(pos, ax=ax, shrink=0.5,location='right')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('TF')
    ax.set_ylabel('TF')

    fig.set_size_inches([14,12])
    plt.tight_layout()
    fig.savefig(outfig)
    plt.close(fig)



