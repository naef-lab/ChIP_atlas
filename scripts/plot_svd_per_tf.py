import numpy as np
import pandas as pd
import h5py
import pickle
import os
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--tf'
        ,required=True
        ,type=str
        ,help="TF")
    parser.add_argument('--infile'
        ,required=True
        ,type=str
        ,help="hdf5 with SVD and pearson's correlation")
    parser.add_argument('--outfig'
        ,required=True
        ,type=str
        ,help="pdf figure")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_argument()
    
    # open hdf5 file
    with h5py.File(args.infile,'r') as hf:
        rho = hf['rho'][:]
        S = hf['S'][:]
        Vh = hf['Vh'][:]

    # experiment hierarchical clustering based on corr matrix
    n = rho.shape[0]
    triu_idx = np.triu_indices(n,1)
    corr_dist = 1 - rho[triu_idx[0],triu_idx[1]]
    linkage = hierarchy.linkage(corr_dist,optimal_ordering=True)

    # Get pcs, loadings and explained variance
    #X_svd = U @ np.diag(S) @ Vh
    Loadings = np.square(Vh)
    explained_var = np.square(S)
    explained_var /= sum(explained_var)

    # plot results (exp. corr, exp. loadings, Variance explained per comp)
    fig, axes = plt.subplots(1,4)

    # plot hierarchichal clustering based on corr
    ax = axes[0]
    R = hierarchy.dendrogram(linkage, ax=ax, above_threshold_color='y',orientation='left',labels=None)
    ax.set_yticklabels([])

    # reorder
    rho = rho[R['leaves'],:]
    rho = rho[:,R['leaves']]
    Loadings = Loadings[:,R['leaves']]

    # plot corr matrix
    ax = axes[1]
    pos = ax.imshow(rho, cmap='RdBu_r', interpolation=None,vmin=-1,vmax=1,origin='lower')
    fig.colorbar(pos, ax=ax, shrink=0.5,location='right')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('Experiments')
    ax.set_ylabel('Experiments')

    # plot experiment loadings
    ax = axes[2]
    pos = ax.imshow(Vh.T, cmap='RdBu_r', interpolation=None,origin='lower')
    fig.colorbar(pos, ax=ax,shrink=0.5, location='right')
    ax.set_yticklabels([])
    ax.set_xlabel('Components')
    ax.set_ylabel('Experiments')

    # plot fraction of explained variance
    ax = axes[3]
    ax.plot(np.arange(explained_var.shape[0]),explained_var,'.-')
    ax.set_yscale('log')
    ax.set_xlabel('Component')
    ax.set_ylabel('Frac of variance')
    
    # save fig
    fig.set_size_inches([4*12,12])
    plt.tight_layout()
    fig.savefig(args.outfig)
    plt.close(fig)




