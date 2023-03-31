import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, single
import matplotlib.pyplot as plt


# input/output files
infile_tf_peaks_jaccard = 'results/Jaccard_peak_pairs.txt'

# get TF peaks pairs jaccard indices
tf_jaccard = pd.read_csv(infile_tf_peaks_jaccard,sep='\t',index_col=0)
tf_jaccard.loc[:,:] = tf_jaccard.values + tf_jaccard.values.T

N_tf = tf_jaccard.shape[0]

D = 1 - tf_jaccard.values[ np.triu_indices(N_tf,1) ]

Z = single(D)

fig = plt.figure()
dn = dendrogram(Z)

fig.set_size_inches([25,10])
plt.tight_layout()
fig.savefig(f'results/fig/hierarchical_clustering/single_func.pdf')
plt.close(fig)

for method in ['single','complete','average','ward']:
    print(method)

    Z = linkage(D, method=method)

    fig = plt.figure()
    dn = dendrogram(Z)

    fig.set_size_inches([25,10])
    plt.tight_layout()
    fig.savefig(f'results/fig/hierarchical_clustering/{method}.pdf')
    plt.close(fig)