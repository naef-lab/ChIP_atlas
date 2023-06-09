import numpy as np
import pandas as pd
import h5py
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--infile'
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

    with h5py.File(args.infile,'r') as hf:
        X = hf['chip_prom_pos_exp'][:]

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
        hf.create_dataset('rho',data=rho)
        hf.create_dataset('U',data=U)
        hf.create_dataset('S',data=S)
        hf.create_dataset('Vh',data=Vh)