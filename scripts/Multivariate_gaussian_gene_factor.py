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


def plot_covariance_and_precision(covs,precs,outfig):

    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.02, right=0.98)

    # plot the covariances
    vmax = max([covs[i][1].max() for i in range(1,len(covs))])
    for i, (name, this_cov) in enumerate(covs):
        vmax = ( this_cov-np.diag(np.diag(this_cov)) ).max()
        plt.subplot(2, len(covs), i + 1)
        plt.imshow(
            this_cov,
            vmin=-vmax,
            vmax=vmax,
            cmap=plt.cm.RdBu_r,
        )
        plt.xticks(())
        plt.yticks(())
        plt.title("%s covariance" % name)

    # plot the precisions
    vmax = max([precs[i][1].max() for i in range(1,len(precs))])
    for i, (name, this_prec) in enumerate(precs):
        vmax = (this_prec-np.diag(np.diag(this_prec)) ).max()
        ax = plt.subplot(2, len(precs), i + 1 + len(precs))
        plt.imshow(
            np.ma.masked_equal(this_prec, 0),
            vmin=-vmax,
            vmax=vmax,
            cmap=plt.cm.RdBu_r,
        )
        plt.xticks(())
        plt.yticks(())
        plt.title("%s precision" % name)
        if hasattr(ax, "set_facecolor"):
            ax.set_facecolor(".7")
        else:
            ax.set_axis_bgcolor(".7")
    
    #plt.tight_layout()
    plt.savefig(outfig)




def plot_tf_correlation_matrix(X,genome):

    rho = np.corrcoef(X.T)
    for method in ['single','complete','average','weighted','centroid','median','ward']:

        # TF hierarchical clustering based on corr matrix
        triu_idx = np.triu_indices(n_tf,1)
        corr_dist = 1 - rho[triu_idx[0],triu_idx[1]]
        linkage = hierarchy.linkage(corr_dist,optimal_ordering=True,method=method)

        # plot results (exp. corr, exp. loadings, Variance explained per comp)
        fig, axes = plt.subplots(1,2)

        ax = axes[0]
        R = hierarchy.dendrogram(linkage, ax=ax, above_threshold_color='y',orientation='left',labels=None)
        ax.set_yticklabels([])

        # reorder corr matrix
        rho = rho[R['leaves'],:]
        rho = rho[:,R['leaves']]

        ax = axes[1]
        pos = ax.imshow(rho, cmap='RdBu_r', interpolation=None,vmin=-1,vmax=1,origin='lower')
        fig.colorbar(pos, ax=ax, shrink=0.5,location='top')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlabel('TF')
        ax.set_ylabel('TF')

        fig.set_size_inches([14,12])
        plt.tight_layout()
        fig.savefig(f'results/fig/{genome}/correlation_tf_Chip_{method}.pdf')
        plt.close(fig)


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

    # box-cox transform
    if False:
        pc = 1e-6
        lam = np.zeros(n_tf)
        for i in range(n_tf):
            X[:,i], lam[i] = boxcox(X[:,i]+pc)

    # regularize cov. matrix
    #pc = 1e-4
    #X = np.log(X+pc)
    X -= X.mean(axis=0)
    #X /= X.std(axis=0)
    #outfile = 'results/TF_mean.npy'
    #np.save(outfile,mu)

    if False:
        plot_tf_correlation_matrix(X,genome)

    emp_cov = np.dot(X.T, X) / n_prom
    emp_prec = linalg.inv(emp_cov)

    lw = LedoitWolf(assume_centered=True).fit(X)
    lw_cov = lw.covariance_
    lw_prec = lw.precision_
    print(lw.shrinkage_)

    oas_cov,oas_shrinkage = oas(X,assume_centered=True)
    oas_prec = linalg.inv(oas_cov)
    print(oas_shrinkage)

    if False:
        model = GraphicalLassoCV(assume_centered=True)
        model.fit(X)
        gl_cov = model.covariance_
        gl_prec = model.precision_

        covs = [
            ("Empirical", emp_cov),
            ("Ledoit-Wolf", lw_cov),
            ("OAS", oas_cov),
            #("GraphicalLassoCV", gl_cov),
            ]
        precs = [
            ("Empirical", emp_prec),
            ("Ledoit-Wolf", lw_prec),
            ("OAS", oas_prec),
            #("GraphicalLasso", gl_prec),
        ]
        outfig = f'results/fig/covariance_and_precision_{N_bin}.pdf'
        plot_covariance_and_precision(covs,precs,outfig)

    my_cov  = lw_cov
    my_prec = lw_prec

    Gene_subset = ['aminoacyl_tRNA_synthetase','ribosomal_protein_genes','core_circadian_clock_genes']
    #Gene_subset = ['core_circadian_clock_genes']

    # 1D conditional multivariate gaussian
    if True:
        for gene_subset in Gene_subset:
            print(gene_subset)
            
            genes = pd.read_csv(f'resources/{gene_subset}.txt',sep='\t')
            
            gene_idx = np.array([]).astype(int)
            for g in genes[genome].values:
                gene_idx = np.concatenate([gene_idx,promoterome[promoterome.gene == g].index.values])

            gene_idx = np.unique(gene_idx)
            n_gene = len(gene_idx)
            rnd_idx = np.random.randint(0,promoterome.shape[0],n_gene)
            
            P={}
            P['cond'] = np.zeros([n_gene,n_tf])
            P['cond_rnd'] = np.zeros([n_gene,n_tf])
            P['marg'] = np.zeros([n_gene,n_tf])
            P['marg_rnd'] = np.zeros([n_gene,n_tf])

            idx_tf = np.arange(n_tf)
            for j in idx_tf:
                if j % int(n_tf/10)==0:
                    print(f'\t{int(10*j/n_tf)/10}')

                idx_1 = np.array([j])
                idx_2 = idx_tf [ idx_tf != idx_1 ]

                S_11 = my_cov[idx_1,idx_1]
                S_12 = my_cov[np.ix_(idx_1,idx_2)]
                S_21 = my_cov[np.ix_(idx_2,idx_1)]
                S_22_inv = linalg.inv( my_cov[np.ix_(idx_2,idx_2)] )
                S_cond = S_11 - S_12 @ S_22_inv @ S_21
                #L_11_inv = 1/my_prec[np.ix_(idx_1,idx_1)]
                #L_12 = my_prec[np.ix_(idx_1,idx_2)]
                #S_cond = L_11_inv
            
                x_1 = X[np.ix_(gene_idx,idx_1)]
                x_2 = X[np.ix_(gene_idx,idx_2)]

                # mu_cond = mu_1 + S_12 @ S_22_inv @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                mu_cond =  S_12 @ S_22_inv @ x_2.T
                # mu_cond = mu_1 - L_11_inv @ L_12 @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                # mu_cond = L_11_inv @ L_12 @ x_2

                for i,x_g1 in enumerate(x_1):
                    P['cond'][i,j] = multivariate_normal.cdf(x_g1, mean=mu_cond[0,i], cov=S_cond)
                    P['marg'][i,j] = multivariate_normal.cdf(x_g1,mean=0,cov=S_11)

                x_rnd_1 = X[np.ix_(rnd_idx,idx_1)]
                x_rnd_2 = X[np.ix_(rnd_idx,idx_2)]
                mu_cond = S_12 @ S_22_inv @ x_rnd_2.T
                #mu_cond = L_11_inv @ L_12 @ x_rnd_2

                for i,x_g1 in enumerate(x_rnd_1):
                    P['cond_rnd'][i,j] = multivariate_normal.cdf(x_g1, mean=mu_cond[0,i], cov=S_cond)
                    P['marg_rnd'][i,j] = multivariate_normal.cdf(x_g1,mean=0,cov=S_11)

            
            fig, axes = plt.subplots(4,1)
            axes = axes.reshape(2*2)
            # gaussian kernel
            GK = {}
            TFs_out = {}
            y = np.linspace(0,1,201)
            y_tick = np.linspace(0,1,5)
            
            for f,p in enumerate(P):
                print(p)
                GK[p] = np.zeros([y.shape[0],n_tf])
                for tf in range(n_tf):
                    GK[p][:,tf] = gaussian_kde(P[p][:,tf])(y)

                GK[p] = GK[p]/GK[p].max(0)
                sort_idx = np.argsort(np.argmax(GK[p],0))
                GK[p] = GK[p][:,sort_idx]

                ax = axes[f]
                ax.imshow(GK[p],origin='lower',interpolation='none',cmap='YlOrBr',aspect='auto',extent=[0-.5,n_tf-.5,0,1])
                ax.set_yticks(y_tick)
                ax.set_yticklabels(y_tick)
                ax.set_xticks([])

                # annotate high tfs
                idx_tf = np.argmax(GK[p],0) > .8*len(y)
                my_tfs = TFs[sort_idx][idx_tf]
                tick_pos = np.arange(n_tf)[idx_tf]
                label_pos = np.linspace(0,n_tf-1,sum(idx_tf))
                y_pos = 1.1
                for i,j in enumerate(np.where(idx_tf)[0]):
                    if sum(idx_tf)>20:
                        ax.text(label_pos[i],y_pos,TFs[sort_idx][j],va='bottom',ha='center',rotation=90)
                    else:
                        ax.text(label_pos[i],y_pos,TFs[sort_idx][j],va='bottom',ha='center')
                    ax.plot([label_pos[i],tick_pos[i]] ,[y_pos,1],'k-',lw=.1)

                # annotate low tfs
                idx_tf = np.argmax(GK[p],0) < .2*len(y)
                my_tfs = TFs[sort_idx][idx_tf]
                tick_pos = np.arange(n_tf)[idx_tf]
                label_pos = np.linspace(0,n_tf-1,sum(idx_tf))
                y_pos = -.1
                for i,j in enumerate(np.where(idx_tf)[0]):
                    if sum(idx_tf)>20:
                        ax.text(label_pos[i],y_pos,TFs[sort_idx][j],va='top',ha='center',rotation=90)
                    else:
                        ax.text(label_pos[i],y_pos,TFs[sort_idx][j],va='top',ha='center')
                    ax.plot([label_pos[i],tick_pos[i]] ,[y_pos,0],'k-',lw=.1)
                if 'rnd' in p:
                    ax.set_ylabel(f'{p} ({n_gene})')
                else:
                    ax.set_ylabel(f'{p} {gene_subset} ({n_gene})')
                plt.rcParams["axes.edgecolor"]=(1,1,1,0)

                if p == 'cond':
                    MY_TFS = my_tfs

            fig.set_size_inches([30,20])
            plt.tight_layout()
            fig.savefig(f'results/fig/{genome}/Conditional_multivariate_marginal_{gene_subset}_{N_bin}.pdf')
            plt.close(fig)


            fig, axes = plt.subplots(2,len(MY_TFS))

            idx_tf = np.arange(n_tf)
            for t,tf in enumerate(MY_TFS):
                idx_1 = np.where(tf==TFs)[0]
                idx_2 = idx_tf [ idx_tf != idx_1 ]

                S_11 = my_cov[idx_1,idx_1]
                S_12 = my_cov[np.ix_(idx_1,idx_2)]
                S_21 = my_cov[np.ix_(idx_2,idx_1)]
                S_22_inv = linalg.inv( my_cov[np.ix_(idx_2,idx_2)] )
                S_cond = S_11 - S_12 @ S_22_inv @ S_21
            
                x_1 = X[np.ix_(gene_idx,idx_1)]
                x_2 = X[np.ix_(gene_idx,idx_2)]

                # mu_cond = mu_1 + S_12 @ S_22_inv @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                mu_cond =  S_12 @ S_22_inv @ x_2.T

                # plot histogram of X[:,tf] and marginal distribution
                ax = axes[0,t]
                ax.hist(X[:,idx_1],50,density=True)
                ax.hist(X[gene_idx,idx_1],50,density=True)

                x = np.linspace( 1.1*X[:,idx_1].min(), 1.1*X[:,idx_1].max(),1000)

                p_x = multivariate_normal.pdf(x,mean=0,cov=S_11)
                ax.plot(x,p_x)
                ax.set_xlim(x[0],x[-1])
                ax.set_title(tf)

                # plot marginal distribution fo each gene
                ax = axes[1,t]

                P_cond = np.zeros([n_gene,len(x)])
                for g,x_g1 in enumerate(x_1):
                    P_cond[g,:] = multivariate_normal.pdf(x, mean=mu_cond[0,g], cov=S_cond)

                ax.imshow(P_cond,origin='lower',extent=[x[0],x[-1],-.5,n_gene-.5],aspect='auto',interpolation='none')
                ax.plot(x_1,np.arange(n_gene),'r.')
                ax.set_yticks(range(n_gene))
                ax.set_yticklabels(promoterome.gene[gene_idx].values)
                ax.set_xlim(x[0],x[-1])
            
            fig.set_size_inches([len(MY_TFS)*10,12])
            plt.tight_layout()
            fig.savefig(f'results/fig/{genome}/Conditional_multivariate_marginal_{gene_subset}_{N_bin}_outlier_TFs.pdf')
            plt.close(fig)




            if False:
                p_val = {}
                sorted_TF = {}
                for i in P:
                    p_val[i] = P[i]
                    
                    sign = np.sign(p_val[i]-.5)
                    p_val[i][sign==1] = 1 - p_val[i][sign==1]
                    p_val[i][p_val[i]<1e-5] = 1e-5
                    p_val[i] = sign * -np.log10(p_val[i])

                    idx = np.argsort(p_val[i].mean(0))
                    p_val[i] = p_val[i][:,idx]
                    #sorted_TFs[i] = TFs[idx]

            if False:
                p = P['cond_rnd'].mean(0)
                dp = P['cond_rnd'].std(0)
                idx = np.argsort(p)
                p = p[idx]
                dp = dp[idx]
                #ax.errorbar(np.arange(n_tf),p,dp,ls='none',label='random subset',lw=1)
                ax.violinplot(P['cond_rnd'])


            if False:
                p = P['cond'].mean(0)
                dp = P['cond'].std(0)
                idx = np.argsort(p)
                p = p[idx]
                dp = dp[idx]
                sorted_TFs = TFs[idx]

                #ax.errorbar(np.arange(n_tf),p,dp,ls='none',label=gene_subset,lw=.5)
                ax.boxplot(P['cond'])
                ax.plot(np.arange(n_tf),p,'k.',ms=1,label=r'$\mu$')
                ax.plot([0,n_tf],[.5,.5],'k-',label='0.5')
            
            #ax.violinplot(p_val['cond'])

            if False:
                idx_txt = np.where( np.abs(p - 0.5) - dp > 0 )[0]
                idx_txt = idx_txt[ np.argsort( - ( np.abs(p[idx_txt] - 0.5) - dp[idx_txt] ) ) ]
                #idx_txt = np.arange(len(p))[len(p)-10:]

                pos_r = n_tf
                pos_l = 0
                for i in idx_txt:
                    if (n_tf - i) < i:
                        x = pos_r
                        pos_r -= 15
                        ha = 'right'
                    else:
                        x = pos_l
                        pos_l += 15
                        ha = 'left'
                    ax.text(x,p[i],sorted_TFs[i],ha=ha,va='center',fontsize=6)

                if gene_subset=='ribosomal_protein_genes':
                    for g in ['Gabpa','Sp1','Yy1']:
                        i = np.where( sorted_TFs == g )[0]
                        ax.plot(i,p[i],'w.',ms=2)
                        ax.text(i,p[i]-3*dp[i],g,ha='center',va='top',fontsize=6,color='blue')

            if False:
                ax.legend(loc='upper left')
                ax.set_xlim([-1,n_tf])
                ax.set_xticklabels([])
                #ax.set_yticklabels([])
                ax.set_xlabel('sorted TF')
                ax.set_ylabel(r'$\int_{-\infty}^{x_1}dx P(x|x_2)$')
                ax.set_title(f'P conditional {gene_subset} ({n_gene})')


                ax = axes[1]

                p = P['marg_rnd'].mean(0)
                dp = P['marg_rnd'].std(0)
                idx = np.argsort(p)
                p = p[idx]
                dp = dp[idx]
                ax.errorbar(np.arange(n_tf),p,dp,ls='none',label='random subset',lw=1)

                p = P['marg'].mean(0)
                dp = P['marg'].std(0)
                idx = np.argsort(p)
                p = p[idx]
                dp = dp[idx]
                sorted_TFs = TFs[idx]

                ax.errorbar(np.arange(n_tf),p,dp,ls='none',label=gene_subset,lw=.5)
                ax.plot(np.arange(n_tf),p,'k.',ms=1,label=r'$\mu$')
                ax.plot([0,n_tf],[.5,.5],'k-',label='0.5')

                idx_txt = np.where( np.abs(p - 0.5) - dp > 0.1 )[0]
                idx_txt = idx_txt[ np.argsort( - ( np.abs(p[idx_txt] - 0.5) - dp[idx_txt] ) ) ]

                pos_r = n_tf
                pos_l = 0
                for i in idx_txt:
                    if (n_tf - i) < i:
                        x = pos_r
                        pos_r -= 15
                        ha = 'right'
                    else:
                        x = pos_l
                        pos_l += 15
                        ha = 'left'
                    ax.text(x,p[i],sorted_TFs[i],ha=ha,va='center',fontsize=6)

                if gene_subset=='ribosomal_protein_genes':
                    for g in ['Gabpa','Sp1','Yy1']:
                        i = np.where( sorted_TFs == g )[0]
                        ax.plot(i,p[i],'w.',ms=2)
                        ax.text(i,p[i]-3*dp[i],g,ha='center',va='top',fontsize=6,color='blue')

                ax.set_xlim([-1,n_tf])
                ax.set_xticklabels([])
                #ax.set_yticklabels([])
                ax.set_xlabel('sorted TF')
                ax.set_ylabel(r'$\int_{-\infty}^{x_1}dx P(x)$')
                ax.set_title(f'P marginal {gene_subset} ({n_gene})')


                fig.set_size_inches([6,8])
                plt.tight_layout()
                fig.savefig(f'results/fig/{genome}/Conditional_multivariate_marginal_{gene_subset}_{N_bin}.pdf')
                plt.close(fig)


    if False:
        Gene_subset = ['mm10_Aminoacyl_tRNA_synthetase','mm10_ribosomal_proteins']
        for gene_subset in Gene_subset:
            print(gene_subset)
            
            genes = pd.read_csv(f'resources/{gene_subset}.txt',sep='\t')
            
            gene_idx = np.array([]).astype(int)
            for g in genes['Gene name'].values:
                gene_idx = np.concatenate([gene_idx,promoterome[promoterome.gene == g].index.values])

            gene_idx = np.unique(gene_idx)
            rnd_idx = np.random.randint(0,promoterome.shape[0],n_gene)


            P = np.zeros([n_gene,n_tf,n_tf])
            P_rnd = np.zeros([n_gene,n_tf,n_tf])
            idx_tf = np.arange(n_tf)

            for j in range(n_tf-1):
                print(f'\t{j/n_tf}')

                for k in range(j+1,n_tf):

                    idx_1 = np.array([j,k])
                    idx_2 = np.array(list(set(idx_tf) - set(idx_1)))
                    
                    S_11 = lw_cov[np.ix_(idx_1,idx_1)]
                    S_12 = lw_cov[np.ix_(idx_1,idx_2)]
                    S_21 = lw_cov[np.ix_(idx_2,idx_1)]
                    S_22_inv = linalg.pinv( lw_cov[np.ix_(idx_2,idx_2)] )

                    S_cond = S_11 - S_12 @ S_22_inv @ S_21
                    
                    
                    for i,g in enumerate( zip(gene_idx,rnd_idx) ): 
                        x_1 = X[g[0],idx_1]
                        x_2 = X[g[0],idx_2]
                        # mu_cond = mu_1 + S_12 @ S_22_inv @ (x_2 - mu_2) with mu_1 = mu_2 = 0 :
                        mu_cond = S_12 @ S_22_inv @ x_2

                        # Mahalanobis distance
                        if True:
                            try:
                                d_M = mahalanobis(x_1,mu_cond,linalg.inv(S_cond))
                            except:
                                d_M = 0
                            P[i,j,k] = chi2.cdf(d_M**2, df=2)
                        else:
                            P[i,j,k] = multivariate_normal.cdf(x_1, mean=mu_cond, cov=S_cond)

                        P[i,k,j] = P[i,j,k]

                        x_rnd_1 = X[g[1],idx_1]
                        x_rnd_2 = X[g[1],idx_2]
                        mu_cond = S_12 @ S_22_inv @ x_rnd_2
                        
                        if True:
                            try:
                                d_M = mahalanobis(x_rnd_1,mu_cond,linalg.inv(S_cond))
                            except:
                                d_M = 0
                            P_rnd[i,j,k] = chi2.cdf(d_M**2, df=2)
                        else:
                            P_rnd[i,j,k] = multivariate_normal.cdf(x_rnd_1, mean=mu_cond, cov=S_cond)
                        P_rnd[i,k,j] = P_rnd[i,j,k]


            fig, axes = plt.subplots(1,2)

            p = P_rnd.mean(0)
            dp = P_rnd.std(0)
            idx = np.argsort(p.max(0))
            p = p[np.ix_(idx,idx)]
            
            ax = axes[0]
            pos = ax.imshow(p, cmap='jet', interpolation=None,vmin=0,vmax=p.max(),origin='lower')
            fig.colorbar(pos, ax=ax, shrink=0.5,location='top',label=r'$\int_0^{d^2} dx \chi^2_k(x)$')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title('random subset')
            ax.set_axis_off()

            p = P.mean(0)
            dp = P.std(0)
            idx = np.argsort(p.max(0))
            p = p[np.ix_(idx,idx)]
            sorted_TFs = TFs[idx]
            
            ax = axes[1]
            pos = ax.imshow(p, cmap='jet', interpolation=None,vmin=0,vmax=p.max(),origin='lower')
            fig.colorbar(pos, ax=ax, shrink=0.5,location='top',label=r'$\int_0^{d^2} dx \chi^2_k(x)$')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(gene_subset)
            ax.set_axis_off()

            fig.set_size_inches([12,6])
            plt.tight_layout()
            fig.savefig(f'results/fig/Conditional_multivariate_2_TFs_{gene_subset}_{N_bin}.pdf')
            plt.close(fig)

            fig, ax = plt.subplots(1,1)
                    
            pos = ax.imshow(p[-10:,-10:], cmap='jet', interpolation=None,vmin=0,vmax=p.max(),origin='lower')
            fig.colorbar(pos, ax=ax, shrink=0.5,location='top',label=r'$\int_0^{d^2} dx \chi^2_k(x)$')
            ax.set_xticks(np.arange(10))
            ax.set_yticks(np.arange(10))
            ax.set_xticklabels(sorted_TFs[-10:])
            ax.set_yticklabels(sorted_TFs[-10:])
            ax.set_title(gene_subset)

            fig.set_size_inches([6,6])
            plt.tight_layout()
            fig.savefig(f'results/fig/Conditional_multivariate_2_TFs_{gene_subset}_{N_bin}_top_10.pdf')
            plt.close(fig)