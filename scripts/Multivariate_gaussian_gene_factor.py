import numpy as np
import pandas as pd
import os
import h5py
from scipy.linalg import pinvh

from scipy.stats import multivariate_normal

if __name__ == '__main__':

    # Get all TFs
    infile='resources/experimentList_mm10_TFs.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()
    TFs[TFs!='Ctcf']

    N_prom = 30114
    N_pos = 100
    N_tf = TFs.shape[0]
    N_pc = 3

    X = np.zeros([N_tf,N_prom*N_pos])

    for t,tf in enumerate(TFs):
        print(np.round(t/N_tf,3))
        
        infiles_svd={'U':f'results/svd/{tf}_U.npy','S':f'results/svd/{tf}_S.npy','Vh':f'results/svd/{tf}_Vh.npy','rho':f'results/corr/{tf}_rho.npy'}
        infile_X=f'results/TF_tensors/{tf}.hdf5'

        if np.all([os.path.exists(f) for f in infiles_svd.values()]):
            U = np.load(infiles_svd['U'])
            S = np.load(infiles_svd['S'])
            if np.sum(U[:,0])<0:
                U[:,0] *= -1
            X[t,:] = U[:,0] * S[0]

        else:
            with h5py.File(infile_X,'r') as hf:
                X[t,:] = hf[tf][:][:,:,0].reshape([N_prom*N_pos])

    # remove TFs with nans everywhere
    idx_out = np.isnan(X).all(axis=1)
    X = X[~idx_out,:]
    TFs = TFs[~idx_out]

    # compute mean, covariance matrix and save
    mu = np.nanmean(X,axis=1)
    outfile = 'results/TF_mean.npy'
    np.save(outfile,mu)
    
    if True:
        # replace nans with 0s
        X[np.isnan(X)] = 0
        C = np.cov(X)
    else:
        C = np.ma.cov(np.ma.masked_invalid(X))
    
    outfile = 'results/TF_covariance_matrix.npy'
    np.save(outfile,C)

    # plot covariance matrix

    # invert covariance
    Sigma = pinvh(C)
    outfile = 'results/TF_inverse_covariance_matrix.npy'
    np.save(outfile,Sigma)

    # my_gene = ...

    # P_x = multivariate_normal(x,mu,Sigma)