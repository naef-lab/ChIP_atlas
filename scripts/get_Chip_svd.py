import numpy as np
import pandas as pd
import h5py
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Compute svd of Chip signal accross experiments')
    parser.add_argument('--infile_chip'
        ,required=True
        ,type=str
        ,help="TF Chip signal")
    parser.add_argument('--infile_peak'
        ,required=True
        ,type=str
        ,help="TF Chip signal")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="hdf5 with SVD and pearson's correlation")
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,help="window size (in kb)")
    parser.add_argument('--bin_size'
        ,default=20
        ,type=int
        ,help="bin size")
    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()

    # get chip signal
    with h5py.File(args.infile_chip,'r') as hf:
        X = hf['chip_prom_pos_exp'][:]
    
    # reshape: concatenate all promoter together
    N_prom, N_pos, N_exp = X.shape
    X = X.reshape([N_prom*N_pos,N_exp])

    # get peaks
    with h5py.File(args.infile_peak,'r') as hf:
        Peaks = hf['peak_prom_pos_exp'][:]
    
    # bin peaks in args.bin_size bins
    Peaks = Peaks.reshape([N_prom,N_pos,args.bin_size,N_exp]).any(axis=2)

    # reshape: concatenate all promoter together
    Peaks = Peaks.reshape([N_prom*N_pos,N_exp])

    # normalize as z score (w.r.t. BG) in log-space
    for i in range(N_exp):

        # get background: signal in non-peak & non-nan regions
        x_bg = X[~Peaks[:,i],i]
        x_bg = np.log( x_bg[~np.isnan(x_bg)] )

        # replace non-nan signal with normalized signal
        idx_data = ~np.isnan(X[:,i])
        idx_nan = ~idx_data
        X[idx_data,i] = (np.log(X[idx_data,i]) - np.mean(x_bg)) / np.std(x_bg)
        X[idx_nan,i] = np.mean(X[idx_data,i])

    # Normalize by size factor
    # X = ( X / np.nansum(X,axis=0,keepdims=True).sum(axis=1,keepdims=True) ) * args.size_factor
    # replace nans with 0s
    # X[np.isnan(X)] = 0

    #idx_svd = np.where((np.isnan(X).sum(axis=1)<=(N_exp-2)))[0] # at least 2 non-nan values
    idx_svd = np.where((Peaks.sum(axis=1)>0))[0] # at least 1 peak

    if len(idx_svd) == 0:
        idx_svd = (X > X.mean(0) + 3*X.std(0)).any(1)

    if N_exp > 1:
        # pearson corr.
        rho = np.corrcoef(X.T)
        # SVD
        U,S,Vh = np.linalg.svd(X[idx_svd,:],full_matrices=False)

        # change sign such that sum(U,0) > 0
        sign = np.sign(U.sum(axis=0,keepdims=True))
        U = sign*U
        Vh = sign.T*Vh

        U_ = np.zeros(X.shape)*np.nan
        U_[idx_svd,:] = U
        U = U_
    
    else:
        rho = np.ones([N_exp,N_exp])
        U = np.zeros([N_prom*N_pos,N_exp])*np.nan
        S = np.ones([N_exp])
        Vh = np.ones([N_exp,N_exp])

        U[idx_svd,:] = X[idx_svd,:]
        # normalize U
        U = U / np.linalg.norm(U,axis=0,keepdims=True)



    # save in output hdf5 file
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset('rho',data=rho)
        hf.create_dataset('U',data=U)
        hf.create_dataset('S',data=S)
        hf.create_dataset('Vh',data=Vh)
    