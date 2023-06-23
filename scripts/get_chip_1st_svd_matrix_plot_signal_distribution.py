import numpy as np
import pandas as pd
import os
import h5py
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # args
    Genome = ['hg38','mm10']
    window_kb = '2'
    for genome in Genome[:1]:
        
        outfile = f'results/{genome}/tensor_TFsvd1_posbin_prom.hdf5'
        outfig = f'results/fig/{genome}/mean_chip_signal_per_promoter_per_tf.pdf'
        
        # get promoterome
        infile_promoter = f"/home/jbreda/Promoterome/results/{genome}/promoterome_pm{window_kb}kb_filtered.bed"
        promoterome = pd.read_csv(infile_promoter ,sep='\t')
        idx_minus_strand = (promoterome.strand=='-').values

        # Get all TFs
        infile=f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab'
        experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3])
        experiment_tf.columns = ['id','antigen']
        TFs = list(experiment_tf.antigen.unique())

        # get dims
        N_prom = promoterome.shape[0]
        N_pos = 100
        N_tf = len(TFs)
        
        # initilize tf x pos x prom tensor
        X = np.zeros([N_tf,N_pos,N_prom])
        for t,tf in enumerate(TFs):
            print(np.round(t/N_tf,3))
            
            infile = f'results/{genome}/svd/{tf}.hdf5'
            with h5py.File(infile,'r') as hf:
                u = hf['U'][:,0]
                s = hf['S'][0]
            
            X[t,:,:] = np.reshape(u*s,[N_prom,N_pos]).T
        
        # flip promoters on minus strans
        X[:,:,idx_minus_strand] = X[:,::-1,idx_minus_strand]
        Chip_signal = X.mean(axis=2).mean(axis=0)

        # get Chip signal cdf
        cdf = np.cumsum(Chip_signal)/sum(Chip_signal)

        # bin X accross position in equal intervals of the cdf
        BIN_IDX = []
        with h5py.File(outfile,'w') as hf:
            for N_bin in range(1,11):
                
                print(N_bin)

                if N_bin==1:
                    X_binned = X.mean(axis=1)
                else:
                    q = np.linspace(0,1,N_bin+1)[:-1]
                    bin_idx  = np.argmin(np.abs( cdf[:,None] - q[None,:] ),axis=0)

                    X_binned = np.zeros([N_tf,N_bin,N_prom])
                    for i in range(N_bin-1):
                        X_binned[:,i,:] = X[:,bin_idx[i]:bin_idx[i+1],:].mean(axis=1)
                    X_binned[:,N_bin-1,:] = X[:,bin_idx[N_bin-1]:,:].mean(axis=1)

                    X_binned = X_binned.reshape([N_tf*N_bin,N_prom])
                    
                    BIN_IDX.append(bin_idx)

                hf.create_dataset(str(N_bin),data=X_binned)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        #x = np.arange(-int(window_kb)*1000,int(window_kb)*1000,20) + 10
        x = np.linspace(-int(window_kb)*1000,int(window_kb)*1000,N_pos+1)
        x = .5*(x[:-1]+x[1:])
        bin_size = 2*int(window_kb)*1000/N_pos
        ax.bar(x=x,height=Chip_signal,width=bin_size)
        ms = np.linspace(8,4,9)
        y = np.linspace(0,Chip_signal.min()/2,10)[1:]
        for n in range(9):
            ax.plot(x[BIN_IDX[n][1:]],np.repeat(y[n],n+1),'o',ms = ms[n],color=plt.cm.tab10(n+1))
        ax.set_xlabel('Pos. rel. to tss')
        ax.set_ylabel('Mean Chip signal per promoter and tf')
        ax.grid('on')
        fig.set_size_inches([8,6])
        plt.tight_layout()
        fig.savefig(outfig)


        