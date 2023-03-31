import numpy as np
import pandas as pd
import os

mu = []
sigma = []
n = []
l_tot = []
for f in os.listdir('results/bedfiles/'):
    try:
        bed = pd.read_csv(f'results/bedfiles/{f}',sep='\t',header=None)
    except:
        print(f)
        continue
    mu.append(bed.iloc[:,3].mean())
    sigma.append(bed.iloc[:,3].std())
    n.append(bed.shape[0])
    l_tot.append(bed.iloc[:,3].sum())

mu = np.array(mu)[None,:]
n = np.array(n)[None,:]
L = 2730855475
L_prom=1282713

p_overlap = (mu.T+mu)/L
nn = n.T*n
mean_intersection = nn*p_overlap
mean_union = n.T+n - mean_intersection
mean_jaccard = mean_intersection / mean_union
