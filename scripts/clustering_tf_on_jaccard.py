import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.cluster.hierarchy import dendrogram, linkage, single
import matplotlib.pyplot as plt


# input/output files
infile_tf_peaks_jaccard = 'results/Jaccard_peak_pairs.txt'

# get TF peaks pairs jaccard indices
tf_jaccard = pd.read_csv(infile_tf_peaks_jaccard,sep='\t',index_col=0)
tf_jaccard.loc[:,:] = tf_jaccard.values + tf_jaccard.values.T

N_tf = tf_jaccard.shape[0]


# get block matrix

# rearange
graph = csr_matrix(tf_jaccard>.3)
idx = reverse_cuthill_mckee(graph)
tf_jaccard = tf_jaccard.iloc[idx,idx]

# plot whole Jaccard index matrix
fig = plt.figure()
ax = fig.add_subplot(111)
pos = ax.imshow(tf_jaccard, cmap='twilight', interpolation=None)

fig.colorbar(pos, ax=ax, shrink=0.5,location='top')
ax.set_xticks(range(tf_jaccard.shape[0]))
ax.set_xticklabels(tf_jaccard.columns,fontsize=6,rotation=90)
ax.set_yticks(range(tf_jaccard.shape[1]))
ax.set_yticklabels(tf_jaccard.index,fontsize=6)

fig.set_size_inches([24,24])
plt.tight_layout()
fig.savefig('results/fig/jaccard_index_tf_peaks.pdf')
plt.close(fig)


if False:
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