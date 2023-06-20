import numpy as np
import pandas as pd
import os
import h5py
from scipy import linalg
from scipy.stats import multivariate_normal, chi2, norm, gaussian_kde, boxcox
from scipy.spatial.distance import mahalanobis
from sklearn.covariance import GraphicalLassoCV, LedoitWolf, oas, empirical_covariance, shrunk_covariance, log_likelihood
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

if __name__ == '__main__':

    genome = 'mm10'
    window = '2'

    # get promoterome
    infile_promoterome=f'/home/jbreda/Promoterome/results/mm10/promoterome_pm{window}kb_filtered.bed'
    promoterome = pd.read_csv(infile_promoterome,sep='\t')

    # Get all TFs
    infile='resources/experimentList_mm10_TFs_only_QC_filtered.tab'
    experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3],index_col=0)
    TFs = experiment_tf.antigen.unique()

    # get binned chip signal per tf per prom
    infile = f'results/{genome}/tensor_TFsvd1_posbin_prom.hdf5'
    N_bin = 1
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
    #X /= X.std(axis=0)
    #outfile = 'results/TF_mean.npy'
    #np.save(outfile,mu)

    emp_cov = np.dot(X.T, X) / n_prom
    emp_prec = linalg.inv(emp_cov)

    lw = LedoitWolf(assume_centered=True).fit(X)
    lw_cov = lw.covariance_
    lw_prec = lw.precision_
    print(lw.shrinkage_)

    oas_cov,oas_shrinkage = oas(X,assume_centered=True)
    oas_prec = linalg.inv(oas_cov)
    print(oas_shrinkage)

    my_cov = lw_cov
    my_prec = lw_prec

    sigma = np.sqrt(np.diag(my_cov))
    my_corr = my_cov/(sigma[:,np.newaxis]  @ sigma[np.newaxis,:])

    diag = np.sqrt(np.diag(my_prec))
    my_par_corr = -my_prec/(diag[:,np.newaxis]  @ diag[np.newaxis,:])



    # Get TF paris annotated in a complex
    Complex = pd.read_csv('results/TF_Complex.tsv',sep='\t',index_col=0)
    
    idx = -np.ones(len(TFs)).astype(int)
    for i,tf in enumerate(TFs):
        idx[i] = np.where(tf==Complex.columns)[0][0]
    Cplx = Complex.values[np.ix_(idx,idx)]
    Cplx = (Cplx>0).astype(int)

    # TODO: Get TF pairs in the same Jaspar cluster
    Cluster_motif_names = pd.read_csv('../Jaspar/resources/interactive_trees/JASPAR_2022_matrix_clustering_vertebrates_CORE_tables/clusters_motif_names.tab',sep='\t',index_col=0,header=None)
    Cluster_motif_names.columns = ['TFs']
    # Get gene list 
    for c in Cluster_motif_names.index:
        Cluster_motif_names.loc[c,'TFs'] = Cluster_motif_names.loc[c,'TFs'].split(',')
    
    # add mouse gene names to lists
    human_mouse = pd.read_csv('../genome/human_to_mouse.txt',sep='\t')
    # remove nan entries
    idx_out = (human_mouse.loc[:,'Gene name'].values.astype('str')=='nan') | (human_mouse.loc[:,'Mouse gene name'].values.astype('str')=='nan')
    human_mouse = human_mouse.loc[~idx_out,:]
    human_mouse.drop_duplicates(inplace=True)

    TF_cluster_dict = {}
    for c in Cluster_motif_names.index:
        for tf in Cluster_motif_names.loc[c,'TFs']:
            TF_cluster_dict[tf] = c

            if tf in human_mouse.loc[:,'Gene name'].values:
                my_tf = human_mouse.loc[human_mouse.loc[:,'Gene name'].values == tf,'Mouse gene name'].values[0]

                if (not my_tf in TF_cluster_dict) and (not my_tf in Cluster_motif_names.loc[c,'TFs']):
                    TF_cluster_dict[my_tf] = c

    TF_clusters = np.zeros(n_tf).astype(int)
    for i,tf in enumerate(TFs):
        if tf in TF_cluster_dict:
            TF_clusters[i] = int(TF_cluster_dict[tf].split('_')[1])

    same_cluster = (TF_clusters[:,None] == TF_clusters[None,:]).astype(int) - np.eye(n_tf).astype(int)
    

    fig,axes = plt.subplots(1,1)

    ax = axes

    ax.plot([-.2,1],[-.2,1],'k',lw=.1)

    x = (my_corr-np.diag(np.diag(my_corr)) )
    x = np.reshape(x,np.prod(x.shape))
    y = my_par_corr-np.diag(np.diag(my_par_corr))
    y = np.reshape(y,np.prod(y.shape))

    col_cplx = np.reshape(Cplx,np.prod(Cplx.shape))
    col_jaspar = np.reshape(same_cluster,np.prod(same_cluster.shape))
    col = 2*col_cplx + col_jaspar
    my_colors = ['grey','red','green','orange']

    idx_perm = np.random.permutation(np.arange(x.shape[0]))[:100]
    #xy = np.vstack([x,y])
    #xy_subset = np.vstack([x[idx_perm],y[idx_perm]])
    #z = np.ones(shape[x])

    for c in np.unique(col):
        idx = np.where(col==c)[0]
        ax.scatter(x[idx],y[idx],c=my_colors[int(c)],s=2,rasterized=True)
    ax.set_xlabel('Correlation')
    ax.set_ylabel('Partial correlation')
    ax.legend(['','none','same motif cluster','known complex','both'])

    [I,J] = np.where(my_par_corr>0.15)
    for i,j in zip(I,J):
        if i<j:
            ax.text(my_corr[i,j],my_par_corr[i,j],f'{TFs[i]}-{TFs[j]}',fontsize=4,ha='center',va='bottom')

    fig.set_size_inches([8,6])
    plt.tight_layout()
    fig.savefig(f'results/fig/{genome}/scatter_covariance_precision_{N_bin}.pdf',dpi=150)
    plt.close(fig)
    

