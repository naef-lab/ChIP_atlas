import numpy as np
import pandas as pd
import h5py
import pickle
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Get all TFs
    infile='resources/experimentList_mm10_TFs.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
    TFs = experiment_tf.loc[:,3].unique()

    N_prom = 30114
    N_pos = 100
    N_tf = TFs.shape[0]
    N_pc = 10

    outfile_1_exp_tf = 'results/TF_with_1_exp_only.txt'

    with open(outfile_1_exp_tf,'w') as f_out:

        for tf in TFs:
            
            infile_X=f'results/TF_tensors/{tf}.hdf5'
            infile_failed_exp=f'results/TF_tensors/{tf}_failed.txt'
            outfiles={'U':f'results/svd/{tf}_U.npy','S':f'results/svd/{tf}_S.npy','Vh':f'results/svd/{tf}_Vh.npy','rho':f'results/corr/{tf}_rho.npy'}
            outfig=f'results/svd/{tf}.npy'

            if np.any([not os.path.exists(f) for f in outfiles.values()]) and (tf != 'Ctcf'):
                print(tf)

                with h5py.File(infile_X,'r') as hf:
                    X = hf[tf][:]

                # remove experiments failing to open
                with open(infile_failed_exp,'r') as f:
                    idx_out = [int(line.strip().split('\t')[0]) for line in f.readlines()]
                if len(idx_out)>0:
                    X = np.delete(X,idx_out,axis=2)

                # replace nans with 0s
                X[np.isnan(X)] = 0
                N_prom, N_pos, N_exp = X.shape

                print(N_exp)
                if N_exp > 1:
                    X = X.reshape([N_prom*N_pos,N_exp])
                    
                    # compute preson corr.
                    rho = np.corrcoef(X.T)
                    np.save(outfiles['rho'],rho)

                    # SVD
                    # mu = np.mean(X,axis=0,keepdims=True)
                    # sigma = np.std(X,axis=0,keepdims=True)
                    # X_norm = (X-mu)/sigma
                    U,S,Vh = np.linalg.svd(X,full_matrices=False)
                    # save X_svd
                    np.save(outfiles['U'],U)
                    np.save(outfiles['S'],S)
                    np.save(outfiles['Vh'],Vh)
                else:
                    f_out.write(f'{tf}\n')

        

