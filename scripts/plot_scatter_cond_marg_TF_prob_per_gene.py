import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

def plot_scatter_Z_score_cond_marg(Z_score,outfig,genome,N_bin,window):

   fig, axes = plt.subplots(1,2,figsize=(2*8,6))

   ax = axes[0]
   ax.scatter(Z_score['cond'].flatten(),Z_score['marg'].flatten(),s=.5,alpha=1,rasterized=True)
   ax.set_xlabel('Z_score conditional')
   ax.set_ylabel('Z_score marginal')
   ax.set_title(f'{genome} - {N_bin} bin - {window}kb window')

   ax = axes[1]
   ax.hist(Z_score['cond'].flatten()-Z_score['marg'].flatten(),bins=100,alpha=.5,density=True)
   ax.set_xscale('symlog')
   ax.set_xlabel('Z_score conditional - Z_score marginal')
   ax.set_ylabel('Density')

   plt.tight_layout()
   fig.savefig(outfig,dpi=300)
   plt.close()



if __name__ == '__main__':
     
   Genome = ['hg38','mm10']
   window = '2'
   N_bin = 1

   my_gene = {'mm10':'Vldlr',
              'hg38':'VLDLR'}

   for genome in Genome[:1]:
      print(genome)

      # get promoterome
      infile_promoterome=f'/home/jbreda/Promoterome/results/{genome}/promoterome_pm{window}kb_filtered.bed'
      promoterome = pd.read_csv(infile_promoterome,sep='\t')

      # Get all TFs
      infile=f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab'
      experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3],index_col=0)
      TFs = experiment_tf.antigen.unique()

      infile_cond = f'results/{genome}/Z_score_cond_{N_bin}bin_{window}kb.npy'
      infile_marg = f'results/{genome}/Z_score_marg_{N_bin}bin_{window}kb.npy'
      Z_score = {}
      Z_score['cond'] = np.load(infile_cond)
      Z_score['marg'] = np.load(infile_marg)

      # plot Z_score values for Vldlr only
      idx = promoterome.gene == my_gene[genome]
      cond = Z_score['cond'][idx,:].mean(axis=0)
      marg = Z_score['marg'][idx,:].mean(axis=0)

      r = np.sqrt(cond**2+marg**2)
      r_gt_1 = r>1

      fig, ax = plt.subplots(1,1,figsize=(8,6))

      ax.plot([-1,2.5],[-1,2.5],'--',color='gray',alpha=.5,lw=.1)
      ax.scatter(cond[~r_gt_1],marg[~r_gt_1],s=.5,alpha=1,rasterized=True)
      ax.scatter(cond[r_gt_1],marg[r_gt_1],s=5,alpha=1,rasterized=True)
      for i in np.where(r_gt_1)[0] :
         ax.text(cond[i],marg[i],TFs[i],fontsize=6)
      ax.set_xlabel('Z_score conditional')
      ax.set_ylabel('Z_score marginal')
      ax.set_title(f'{genome} - {N_bin} bin - {window}kb window - {my_gene[genome]}')
      if genome == 'hg38':
         plt.axis([-2, 4, -1, 4])

      plt.tight_layout()
      outfig = f'results/fig/{genome}/scatter_Z_score_cond_marg_{N_bin}bin_{window}kb_{my_gene[genome]}.pdf'
      fig.savefig(outfig,dpi=300)
      plt.close()

      # Plot scatter
      if False:
         outfig = f'results/fig/{genome}/scatter_Z_score_cond_marg_{N_bin}bin_{window}kb.pdf'
         plot_scatter_Z_score_cond_marg(Z_score,outfig,genome,N_bin,window)

         # get gene-gene correlation in Zscore diff space
         D = Z_score['cond']-Z_score['marg']
         rho = np.corrcoef(D)
            
         # experiment hierarchical clustering based on corr matrix
         n = rho.shape[0]
         triu_idx = np.triu_indices(n,1)
         corr_dist = 1 - rho[triu_idx[0],triu_idx[1]]
         #linkage = hierarchy.linkage(corr_dist,optimal_ordering=True)
         linkage = hierarchy.ward(corr_dist)

         # Cluster genes based on linkage
         K = hierarchy.fcluster(linkage, t=.1, criterion='distance')
         K_uniq, K_count = np.unique(K,return_counts=True)
         #K = hierarchy.ward(corr_dist)

         my_k = np.where(K_count>100)[0]+1

         for k in my_k:
            if sum(K==k) > 100:
               print(k)
               print(np.sort(promoterome.loc[K==k,'gene'].unique()))

         # plot correlation matrix
         outfig = f'results/fig/{genome}/Correlation_Z_score_cond_marg_{N_bin}bin_{window}kb.pdf'
         fig, axes = plt.subplots(1,2,figsize=(2*8,6))

         ax = axes[0]
         R = hierarchy.dendrogram(linkage, ax=ax, above_threshold_color='y',orientation='left',labels=None)
         ax.set_yticklabels([])

         # reorder
         rho = rho[R['leaves'],:]
         rho = rho[:,R['leaves']]
         
         # plot corr matrix
         ax = axes[1]
         pos = ax.imshow(rho, cmap='RdBu_r', interpolation=None,vmin=-1,vmax=1,origin='lower')
         fig.colorbar(pos, ax=ax, shrink=0.5,location='right')
         ax.set_xticklabels([])
         ax.set_yticklabels([])
         ax.set_xlabel('Gene')
         ax.set_ylabel('Gene')
               
         plt.tight_layout()
         fig.savefig(outfig,dpi=300)
         plt.close()


      

