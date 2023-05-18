import numpy as np
import pandas as pd
import os
import h5py
from scipy.linalg import pinvh
from scipy.stats import multivariate_normal

if __name__ == '__main__':

    # Get all TFs
    infile='resources/experimentList_mm10_TFs_only_QC_filtered.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()

    N_bin = 10

    infile = 'results/tensor_TFsvd1_posbin_prom.hdf5'

    with h5py.File(infile,'r') as hf:
        X = hf[str(N_bin)][:]

    # compute mean, covariance matrix and save
    #mu = np.nanmean(X,axis=1)
    #outfile = 'results/TF_mean.npy'
    #np.save(outfile,mu)

    C = np.cov(X)
    
    #outfile = 'results/TF_covariance_matrix.npy'
    #np.save(outfile,C)

    # plot covariance matrix

    # invert covariance
    Sigma = pinvh(C)
    #outfile = 'results/TF_inverse_covariance_matrix.npy'
    #np.save(outfile,Sigma)

    # my_gene = ...

    # P_x = multivariate_normal(x,mu,Sigma)