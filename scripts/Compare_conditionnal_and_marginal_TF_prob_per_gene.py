import numpy as np
import pandas as pd
import h5py
from scipy import linalg
from sklearn.covariance import GraphicalLassoCV, LedoitWolf, oas, empirical_covariance, shrunk_covariance, log_likelihood


if __name__ == '__main__':

    Genome = ['hg38','mm10']
    window = '2'
    N_bin = 1

    for genome in Genome:
        print(genome)

        # get promoterome
        infile_promoterome=f'/home/jbreda/Promoterome/results/{genome}/promoterome_pm{window}kb_filtered.bed'
        promoterome = pd.read_csv(infile_promoterome,sep='\t')

        # Get all TFs
        infile=f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab'
        experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3],index_col=0)
        TFs = experiment_tf.antigen.unique()

        # get binned chip signal per tf per prom
        infile = f'results/{genome}/tensor_TFsvd1_posbin_prom.hdf5'
        with h5py.File(infile,'r') as hf:
            X = hf[str(N_bin)][:]
        # row: samples (prom_pos)
        # cols: feature (tf)
        X = X.T
        [n_prom, n_tf] = X.shape

        #rename TFs based on binning
        if N_bin > 1:
            TFs = np.array( [f'{tf}_{p}' for tf in TFs for p in range(N_bin)] )

        # regularize cov. matrix
        #pc = 1e-4
        #X = np.log(X+pc)
        X -= X.mean(axis=0)

        #emp_cov = np.dot(X.T, X) / n_prom
        #emp_prec = linalg.inv(emp_cov)

        lw = LedoitWolf(assume_centered=True).fit(X)
        lw_cov = lw.covariance_
        lw_prec = lw.precision_
        print(lw.shrinkage_)

        #oas_cov,oas_shrinkage = oas(X,assume_centered=True)
        #oas_prec = linalg.inv(oas_cov)
        #print(oas_shrinkage)

        my_cov  = lw_cov
        my_prec = lw_prec

        Z_score = {'cond':np.zeros([n_prom,n_tf]),
                'marg':np.zeros([n_prom,n_tf])}
        
        for g in range(n_prom):
            print(f'\t{g}/{n_prom}')
            
            for t in range(n_tf):

                idx_1 = np.array([t])
                idx_2 = np.hstack([np.arange(t),np.arange(t+1,n_tf)])

                S_11 = my_cov[np.ix_(idx_1,idx_1)]
                #S_12 = my_cov[np.ix_(idx_1,idx_2)]
                #S_21 = my_cov[np.ix_(idx_2,idx_1)]
                #S_22_inv = linalg.inv( my_cov[np.ix_(idx_2,idx_2)] )
                #S_cond = S_11 - S_12 @ S_22_inv @ S_21
                L_11_inv = 1/my_prec[np.ix_(idx_1,idx_1)]
                L_12 = my_prec[np.ix_(idx_1,idx_2)]
                S_cond = L_11_inv
            
                x_1 = X[np.ix_([g],idx_1)]
                x_2 = X[np.ix_([g],idx_2)]

                # mu_cond = mu_1 + S_12 @ S_22_inv @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                # mu_cond =  S_12 @ S_22_inv @ x_2.T
                # mu_cond = mu_1 - L_11_inv @ L_12 @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                mu_cond =  - L_11_inv @ L_12 @ x_2.T

                Z_score['cond'][g,t] = (x_1 - mu_cond)/np.sqrt(S_cond)
                Z_score['marg'][g,t] = x_1/np.sqrt(S_11)

        outfile = f'results/{genome}/Z_score_cond_{N_bin}bin_{window}kb.npy'
        np.save(outfile,Z_score['cond'])

        outfile = f'results/{genome}/Z_score_marg_{N_bin}bin_{window}kb.npy'
        np.save(outfile,Z_score['marg'])