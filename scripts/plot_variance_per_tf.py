import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
import sys
sys.path.insert(0, '/home/jbreda/rand_cmap')
from rand_cmap import rand_cmap
sys.path.insert(0, '/home/jbreda/Mixture_models')
from mixture_models import normal_mixture_EM
import pyBigWig


def get_peak_size_per_experiment(exp_id,genome):
    CHR = [f'chr{c+1}' for c in range(19)] + ['chrX','chrY']

    file_in = f'resources/tracks/{genome}/{exp_id}.05.bb'
    with pyBigWig.open(file_in,'r') as bb:
        L = bb.chroms()
        peak_len = []
        peak_m10log10_q = []
        for chr in CHR:
            if chr in L:
                peaks = bb.entries(chr,0,L[chr])
                peak_len.extend( [p[1] - p[0] for p in peaks] )
                peak_m10log10_q.extend( [p[2] for p in peaks] )

    peak_len = np.array(peak_len)
    peak_m10log10_q = np.array(peak_m10log10_q).astype(float)

    return peak_len.mean(), peak_len.std(), peak_m10log10_q.mean(), peak_m10log10_q.std()


def plot_variance_per_tf(genome,N_bin):

    # variance per tf (1st svd binned in promoter)
    infile = f'results/{genome}/tensor_TFsvd1_posbin_prom.hdf5'

    with h5py.File(infile,'r') as hf:
        X = hf[N_bin][:]
    # row: sample (prom)
    # cols: feature (tf)
    X = X.T

    v = X.var(axis=0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(np.log10(v),bins=30,density=True)
    #ax.set_xscale('log',base=2)
    #ax.set_yscale('log')
    ax.set_xlabel(r'$log_{10}$ variance per accross promoters and positions')
    ax.set_ylabel('density')

    #N_exp_per_antigen.plot(kind='bar',ax=ax,width=1)
    fig.set_size_inches([6,3])
    plt.tight_layout()
    fig.savefig(f'results/fig/hist_tf_var_{genome}_{N_bin}bin_per_prom.pdf')

def plot_exp_var_QC(chip,y_cols):

    x_col = '$log_{10} \sigma^2$'
    logLik, rho, mu1, mu2, sigma2, sigma1, ppG1, ppG2 = normal_mixture_EM(chip[x_col],0.5,-2,1,1,1)
    
    fig = plt.figure()
    cols = 4
    rows = 4
    f = 0

    f+=1
    ax = fig.add_subplot(rows,cols,f)
    #ax.hist(chip[x_col].values,bins=100,density=True)
    x_ = np.linspace( chip[x_col].min(), chip[x_col].max(),1000)
    p_ = rho*norm.pdf(x=x_,loc=mu1,scale=sigma1) + (1-rho)*norm.pdf(x=x_,loc=mu2,scale=sigma2)
    ax.plot(x_,p_)
    sep_by_end = [ chip.loc[chip.paired_single_end==i,x_col] for i in chip.paired_single_end.unique() ]
    ax.hist(sep_by_end,bins=50,density=True,rwidth=1)
    ax.legend(['bimodal fit whole dataset','unknown','single end','paired end'])
    ax.set_xlabel(x_col)
    ax.set_ylabel('density')

    f+=1
    ax = fig.add_subplot(rows,cols,f)

    frac_high_mode = np.zeros([len(TFs),2])
    for i,tf in enumerate(TFs):
        frac_high_mode[i,0] = np.sum(ppG2.loc[chip.antigen==tf]>.5)/sum(chip.antigen==tf)
        frac_high_mode[i,1] = sum(chip.antigen==tf)

    ax.plot([0,frac_high_mode[:,1].max()],2*[(1-rho)],'r')

    ax.scatter(frac_high_mode[:,1],frac_high_mode[:,0])
    ax.set_xlabel('nr of exp')
    ax.set_ylabel('frac. in high var mode')
    
    #colors = [np.where(tf==TFs)[0][0] for tf in chip.antigen ]
    #colors = [np.where(i==chip.celltype_class.unique())[0][0] for i in chip.celltype_class ]
    #id_3 = [i[:3] for i in chip.index]
    #colors = [np.where(i==np.unique(id_3))[0][0] for i in id_3 ]
    #cmap = rand_cmap(len(np.unique(colors)), type='bright')
    #cmap = plt.cm.Set1(range(3))
    cmap = ['lightgrey','green','red']
    colors = [cmap[i] for i in chip.paired_single_end]

    chip_sort = chip.sort_values('paired_single_end')

    X = np.zeros([chip.shape[0],len(y_cols)])
    for j,y in enumerate(y_cols):
        f+=1
        ax = fig.add_subplot(rows,cols,f)
        chip_sort.plot.scatter(x_col,y,s=5,c=colors,ax=ax,alpha=.5)
        #ax.errorbar(x=mean_chip[x_col,'mean'],y=mean_chip[y,'mean'],xerr=mean_chip[x_col,'std'],yerr=mean_chip[y,'std'],ls='none')

        X[:,j] = chip.loc[:,y].values

    rho = np.corrcoef( X.T ) - np.eye(X.shape[1])
    vmax = np.abs(rho).max()

    f+=2
    ax = fig.add_subplot(rows,cols,f)
    im = ax.imshow(rho,vmin=-vmax,vmax=vmax,cmap=plt.cm.RdBu_r)
    plt.colorbar(im,fraction=0.046, pad=0.04)
    ax.set_xticks(range(rho.shape[0]))
    ax.set_xticklabels(y_cols,rotation=30)
    ax.set_yticks(range(rho.shape[0]))
    ax.set_yticklabels(y_cols)



    fig.set_size_inches([cols*8,rows*6])
    plt.tight_layout()
    fig.savefig(f'results/fig/{genome}/scatter_exp_var_QC.pdf')

    plt.close('all')


if __name__ == '__main__':


    genome = 'mm10'
    N_bin = '1'
    plot_variance_per_tf(genome,N_bin)

    # get chip experiment table
    infile=f"resources/experimentList_{genome}_TFs_only_QC_filtered.tab"
    chip = pd.read_csv(infile,sep='\t',usecols=[0,1,2,3,4,5,6,7,8])
    chip.columns = ['id','genome','antigen_class','antigen','celltype_class','celltype','celltype_description','QC','title']
    TFs = np.sort(chip.antigen.unique())

    # get promoterome
    infile_promoterome='/home/jbreda/Promoterome/results/mm10/promoterome_pm1kb_filtered.bed'
    promoterome = pd.read_csv(infile_promoterome,sep='\t')
    idx_prom = {'.': (promoterome.black_listed=='.').values,
                'high_signal': (promoterome.black_listed=='High Signal Region').values,
                'low_mappability': (promoterome.black_listed=='Low Mappability').values }

    var_file = f'results/{genome}_var_per_experiment.txt'

    if os.path.exists(var_file):
        Var = pd.read_csv(var_file,sep='\t',index_col=0)
    else:
        IDs = []
        Var = np.zeros([0,3])
        for tf in TFs:
            print(tf)
            
            # get tf ids
            ids = list(chip.loc[chip.antigen==tf,'id'])
            
            # get chip tracks and compute variance
            infile = f'results/{genome}/TF_tensors/{tf}.hdf5'
            with h5py.File(infile,'r') as hf:
                X = hf['chip_prom_pos_exp'][:]
                if 'failed_experiment' in hf.keys():
                    ids.remove( hf['failed_experiment'] )

            [N_prom,N_pos,N_exp] = X.shape
            X = X.reshape([N_prom*N_pos,N_exp])

            IDs.extend(ids)
                      
            #nanvar = np.nanvar(X, axis=0)
            idx_nan = np.isnan(X)
            f_nan = idx_nan.sum(axis=0)/(N_prom*N_pos)
            X[idx_nan] = 0
            var = np.var(X, axis=0)
            mu = np.mean(X, axis=0)

            Var = np.concatenate([Var,np.stack([var,mu,f_nan],axis=1)],axis=0)

        # make table and save
        Var = pd.DataFrame(data=Var,index=IDs,columns=['Var','mu','frac. nan'])
        Var.to_csv(var_file,sep='\t')
    
    # merge chip and log10 var
    chip.set_index('id',drop=True,inplace=True)
    Var.loc[:,['Var','mu']] = np.log10(Var.loc[:,['Var','mu']])
    Var.columns = [r'$log_{10} \sigma^2$',r'$log_{10} \mu$','frac. nan']
    chip = chip.join(Var,how='inner')

    # compare exp variance with QC values
    QC = pd.DataFrame([[float(n) for n in qc.split(',')] for qc in chip.QC],columns=['n_reads','f_mapped','f_duplicates','n_peaks'],index=chip.index)
    QC.loc[:,'f_mapped'] /= 100
    QC.loc[:,'f_duplicates'] /= 100
    QC['peaks_per_unique_mapped_reads'] = QC.n_peaks/(QC.n_reads * QC.f_mapped * (1-QC.f_duplicates))
    QC.loc[(QC.n_peaks==0) | (QC.f_mapped==0),'peaks_per_unique_mapped_reads'] = 0
    # take log 10
    QC.loc[:,'n_reads'] = np.log10( QC.loc[:,'n_reads'] )
    QC.loc[:,'n_peaks'] = np.log10( QC.loc[:,'n_peaks'] + 1 )
    QC.loc[:,'peaks_per_unique_mapped_reads'] = np.log10( QC.loc[:,'peaks_per_unique_mapped_reads'] + 1e-6 )
    QC.columns = [r'$log_{10}$ nr. of reads', 'frac. mapped', 'frac. duplicates', r'$log_{10}$ nr. of peaks + 1',r'$log_{10}$ peaks per read']

    # add to chip table
    chip = pd.concat([chip,QC],axis=1)

    # compute mean and std of float columns
    my_cols = chip.columns[chip.dtypes=='float64']
    mean_chip = chip.groupby('antigen')[my_cols].agg(['mean','std'])

    # get paired-/ single-end experiments
    paired_end = ['paired' in tit for tit in chip.title]
    single_end = [('sequencing' in tit) & ('paired' not in tit) for tit in chip.title]
    paired_single_end = 2*np.array(paired_end).astype(int) + np.array(single_end).astype(int)
    chip.loc[:,'paired_single_end'] = paired_single_end

    # get peak len stats (mu,sigma) and -10log10(q-value) (mu,sigma) per experiment
    peak_stat_file = f'results/{genome}_peak_stats.txt'
    if os.path.exists(peak_stat_file):
        Peak_stats = pd.read_csv(peak_stat_file,sep='\t',index_col=0)
    else:
        print('get peak stats')    
        Peak_stats = np.zeros([chip.shape[0],4])
        for i,exp_id in enumerate(chip.index):
            if (i % chip.shape[0]/100)==0:
                print(i/chip.shape[0])
            Peak_stats[i,:] = get_peak_size_per_experiment(exp_id,genome)
        Peak_stats = pd.DataFrame(Peak_stats,columns=['peak len mean','peak len std','peak -10log10q mean','peak -10log10q std'],index=chip.index)
        Peak_stats.to_csv(peak_stat_file ,sep='\t')
        print('done')
    chip = pd.concat([chip,Peak_stats],axis=1)

    # plot variance per experiment vs QC
    y_cols = list(Var.columns[1:]) + list(QC.columns) + list(Peak_stats.columns)
    plot_exp_var_QC(chip,y_cols)
    
    # get experiments of chip paired-end only, separate by black-listed promoters
    chip_paired = chip.loc[chip.paired_single_end==2,:]
    N_pos = 100
    X_paired = {}
    for idx in idx_prom:
        n_prom = sum(idx_prom[idx])
        X_paired[idx] = np.zeros([chip_paired.shape[0],n_prom*N_pos])
    prom_pos_gt_N = {}
    N = 10
    for c,id in enumerate(chip_paired.index):
        print(c/chip_paired.shape[0])

        # get tf
        tf = chip_paired.loc[id,'antigen']

        # get chip tracks of id
        infile = f'results/{genome}/TF_tensors/{tf}.hdf5'
        with h5py.File(infile,'r') as hf:
            tf_idx = np.where( id == hf['experiment_id'][:].astype(str) )[0][0]
            X = hf['chip_prom_pos_exp'][:,:,tf_idx]

        [Prom,Pos] = np.where(X>N)
        for prom,pos in zip(Prom,Pos):
            if prom not in prom_pos_gt_N:
                prom_pos_gt_N[prom] = []
            prom_pos_gt_N[prom].append(pos)

        for idx in idx_prom:
            [n_prom,n_pos] = X[idx_prom[idx],:].shape
            X_paired[idx][c,:] = X[idx_prom[idx],:].reshape([n_prom*n_pos])

    # get unique prom and pos with extreme value
    for prom in prom_pos_gt_N:
        prom_pos_gt_N[prom] = np.unique(prom_pos_gt_N[prom])

    # get intervals
    for prom in prom_pos_gt_N:
        intervals = []
        for p in prom_pos_gt_N[prom]:
            if len(intervals)==0:
                intervals.append( [p,p] )
            else:
                if p == intervals[-1][1] + 1:
                    intervals[-1][1] = p
                else:
                    intervals[-1] = tuple(intervals[-1])
                    intervals.append( [p,p] )
        intervals[-1] = tuple(intervals[-1])
        prom_pos_gt_N[prom] = intervals

    # add promoter names to prom_pos_gt_N
    for prom in prom_pos_gt_N:
        prom_pos_gt_N[prom].insert(0, promoterome.loc[prom,'id']+'_'+promoterome.loc[prom,'gene']+'_'+promoterome.loc[prom,'black_listed'] )


    # get experiments of chip single-end only, separate by black-listed promoters
    chip_single = chip.loc[chip.paired_single_end==1,:]
    N_prom = 30102
    N_pos = 100
    X_single = {}
    for idx in idx_prom:
        n_prom = sum(idx_prom[idx])
        X_single[idx] = np.zeros([chip_single.shape[0],n_prom*N_pos])

    for c,id in enumerate(chip_single.index):
        print(c/chip_single.shape[0])

        # get tf
        tf = chip_single.loc[id,'antigen']

        # get chip tracks of id
        infile = f'results/{genome}/TF_tensors/{tf}.hdf5'
        with h5py.File(infile,'r') as hf:
            tf_idx = np.where( id == hf['experiment_id'][:].astype(str) )[0][0]
            X = hf['chip_prom_pos_exp'][:,:,tf_idx]

        for idx in idx_prom:
            [n_prom,n_pos] = X[idx_prom[idx],:].shape
            X_single[idx][c,:] = X[idx_prom[idx],:].reshape([n_prom*n_pos])

    # replace nans by 0
    for idx in idx_prom:
        idx_nan = np.isnan(X_paired[idx])
        X_paired[idx][idx_nan] = 0
        idx_nan = np.isnan(X_single[idx])
        X_single[idx][idx_nan] = 0

    fig = plt.figure()
    cols = 4
    rows = 4
    f = 0
    leg = [f'{k} ({idx_prom[k].sum()})' for k in idx_prom.keys() ]
    for i in range(8):

        to_plot = []
        tit = '|'
        for idx in idx_prom:
            to_plot.append( np.log10(X_single[idx][i,:]+1) )
            tit += f' {idx} var: {np.round(np.log10( np.var(X_single[idx][i,:])) ,1)} |'

        f+=1
        ax = fig.add_subplot(rows,cols,f)
        ax.hist(to_plot,bins=50,density=True)
        ax.set_yscale('log')
        ax.set_xlabel('single-end  '+chip_single.index[i])
        ax.set_title(tit)
        ax.legend(leg)

    for i in range(8):

        to_plot = []
        tit = '|'
        for idx in idx_prom:
            to_plot.append( np.log10(X_paired[idx][i,:]+1) )
            tit += f' {idx} var: {np.round(np.log10( np.var(X_paired[idx][i,:])),1)} |'
        
        f+=1
        ax = fig.add_subplot(rows,cols,f)
        ax.hist(to_plot,bins=50,density=True)
        ax.set_yscale('log')
        ax.set_xlabel('paired-end  '+chip_paired.index[i])
        ax.set_title(tit)
        ax.legend(leg)

    fig.set_size_inches([cols*6,rows*4.5])
    plt.tight_layout()
    fig.savefig(f'results/fig/{genome}/hist_bw_paired_single.pdf')

    plt.close('all')

