import numpy as np
import pandas as pd
import h5py
import pickle
import os
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--tf'
        ,required=True
        ,type=str
        ,help="TF")
    parser.add_argument('--infile_tensor'
        ,required=True
        ,type=str
        ,help="TF tensor")
    parser.add_argument('--infile_failed'
        ,required=True
        ,type=str
        ,help="TF tensor")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="hdf5 with SVD and pearson's correlation")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_argument()

    with h5py.File(args.infile_tensor,'r') as hf:
        X = hf[args.tf][:]

    # remove experiments failing to open
    with open(args.infile_failed,'r') as f:
        idx_out = [int(line.strip().split('\t')[0]) for line in f.readlines()]
    if len(idx_out)>0:
        X = np.delete(X,idx_out,axis=2)

    # replace nans with 0s
    X[np.isnan(X)] = 0
    N_prom, N_pos, N_exp = X.shape

    # reshape: concatenate all promoter together
    X = X.reshape([N_prom*N_pos,N_exp])

    # compute pearson corr.
    rho = np.corrcoef(X.T)
    # SVD
    U,S,Vh = np.linalg.svd(X,full_matrices=False)

    # change sign such that sum(U,0) > 0
    sign = np.sign(U.sum(axis=0,keepdims=True))
    U = sign*U
    Vh = sign.T*Vh

    # open hdf5 file
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset('rho',data=rho,compression="gzip")
        hf.create_dataset('U',data=U,compression="gzip")
        hf.create_dataset('S',data=S,compression="gzip")
        hf.create_dataset('Vh',data=Vh,compression="gzip")