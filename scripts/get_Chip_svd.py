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
    
    # bin peaks in args.bin_size bins and concatenate all promoter together
    Peaks = Peaks.reshape([N_prom,N_pos,args.bin_size,N_exp]).any(axis=2)
    Peaks = Peaks.reshape([N_prom*N_pos,N_exp])

    # normalize as z score (w.r.t. BG) in log-space
    #for i in range(N_exp):
    #
    #    # get background: signal in non-peak & non-nan regions
    #    x_bg = X[~Peaks[:,i],i]
    #    x_bg = np.log( x_bg[~np.isnan(x_bg)] )
    #
    #    # replace non-nan signal with normalized signal
    #    idx_data = ~np.isnan(X[:,i])
    #    idx_nan = ~idx_data
    #    X[idx_data,i] = (np.log(X[idx_data,i]) - np.mean(x_bg)) / np.std(x_bg)
    #    X[idx_nan,i] = np.mean(X[idx_data,i])

    # Normalize by size factor
    # X = ( X / np.nansum(X,axis=0,keepdims=True).sum(axis=1,keepdims=True) ) * args.size_factor

    # replace nans with 0s
    X[np.isnan(X)] = 0

    # Do svd in regrions with at least 1 peak  
    idx_svd = np.where((Peaks.sum(axis=1)>0))[0]

    # if no peak, use regions with high signal (3 std above mean)
    if len(idx_svd) == 0:
        idx_svd = (X > X.mean(0) + 3*X.std(0)).any(1)

    if N_exp > 1:
        # pearson corr.
        rho = np.corrcoef(X.T)

        # Compute SVD on peak regions only
        U_peak,S,Vh = np.linalg.svd(X[idx_svd,:],full_matrices=False)

        # change components sign such that sum(U,0) > 0
        sign = np.sign(U_peak.sum(axis=0,keepdims=True))
        U_peak = sign*U_peak
        Vh = sign.T*Vh

        # Get U for all regions
        U = X @ Vh.T @ np.diag(1/S)
    
    else:
        rho = np.ones([N_exp,N_exp])
        U = X
        S = np.ones([N_exp])
        Vh = np.ones([N_exp,N_exp])

    # normalize U
    U = U / np.sqrt((U**2).sum(axis=0,keepdims=True))

    # save in output hdf5 file
    with h5py.File(args.outfile,'w') as hf:
        hf.create_dataset('rho',data=rho)
        hf.create_dataset('U',data=U)
        hf.create_dataset('S',data=S)
        hf.create_dataset('Vh',data=Vh)
    