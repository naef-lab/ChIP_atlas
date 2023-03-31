import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
import matplotlib.pyplot as plt

# input/output files
infile_complex = 'results/TF_Complex.tsv'
infile_tf_peaks_jaccard = 'results/Jaccard_peak_pairs.txt'

# get Complex Portal interaction matrix
tf_complex = pd.read_csv(infile_complex,sep='\t',index_col=0)

# get TF peaks pairs jaccard indices
tf_jaccard = pd.read_csv(infile_tf_peaks_jaccard,sep='\t',index_col=0)
tf_jaccard.loc[:,:] = tf_jaccard.values + tf_jaccard.values.T

# remove 0 summing columns
tf_complex = tf_complex.loc[tf_complex.sum()>0,tf_complex.sum()>0]

# rearange
graph = csr_matrix(tf_complex)
idx = reverse_cuthill_mckee(graph)
tf_complex = tf_complex.iloc[idx,idx]

# find blocks
blocks_border = []
n = tf_complex.shape[0]
X = tf_complex.values
for i in range(1,n-1):
    if X[:i,i:].sum()==0:
        blocks_border.append(i)

# plot complex matrix
fig = plt.figure()

ax = fig.add_subplot(121)
for i in blocks_border:
    ax.plot([-.5,n],[i-.5,i-.5],'k-',lw=.5,alpha=0.1)
    ax.plot([i-.5,i-.5],[-.5,n],'k-',lw=.5,alpha=0.1)
pos = ax.imshow(np.log10(tf_complex), cmap='jet', interpolation=None)

# write tf names next to first column
for i in range(X.shape[0]):
    last_zero_idx = np.nonzero(X[:,i])[0][0] - 1
    x = min(i-1,last_zero_idx)
    ax.text(x,i,tf_complex.index[i],ha='right',va='center',fontsize=6)
    ax.text(i,x,tf_complex.index[i],ha='center',va='baseline',fontsize=6,rotation=90)

#fig.colorbar(pos, ax=ax, shrink=0.5,location='top')
#ax.set_xticks(range(tf_complex.shape[0]))
#ax.set_xticks([])
#ax.set_xticklabels(tf_complex_complex.columns,fontsize=6,rotation=90)
#ax.set_yticks([])
#ax.set_yticklabels(tf_complex.columns,fontsize=6)
ax.axis('off')

# plot Jaccard index on same TFs
ax = fig.add_subplot(122)
pos = ax.imshow(tf_jaccard.loc[tf_complex.index,tf_complex.columns], cmap='jet', interpolation=None)

fig.colorbar(pos, ax=ax, shrink=0.5,location='top')
ax.set_xticks(range(tf_complex.shape[0]))
ax.set_xticklabels(tf_complex.columns,fontsize=6,rotation=90)
ax.set_yticks(range(tf_complex.shape[1]))
ax.set_yticklabels(tf_complex.index,fontsize=6)

fig.set_size_inches([24,12])
plt.tight_layout()
fig.savefig('results/fig/complex_interaction_matrix_jaccard.pdf')
plt.close(fig)



# plot hist of jaccard index that are in complexes agains others
I,J = np.where( tf_complex>0 )
jac_complex = [tf_jaccard.loc[tf_complex.index,tf_complex.columns].iat[i,j] for i, j in zip(I,J)]

I,J = np.where( tf_complex==0 )
jac_not_complex = [tf_jaccard.loc[tf_complex.index,tf_complex.columns].iat[i,j] for i, j in zip(I,J)]

fig = plt.figure()
ax = fig.add_subplot(111)

[h,bins] = np.histogram(jac_complex,100)
x = .5*(bins[:-1]+bins[1:])
ax.plot(x,h)
[h,bins] = np.histogram(jac_not_complex,100)
x = .5*(bins[:-1]+bins[1:])
ax.plot(x,h)
ax.set_yscale('log')
ax.legend(['anotated complex','no complex'])
ax.set_xlabel('Jaccard index')
ax.set_ylabel('Nr. of tf pairs')

fig.set_size_inches([8,6])
plt.tight_layout()
fig.savefig('results/fig/histogram_jaccard_in_out_complex.pdf')
plt.close(fig)



