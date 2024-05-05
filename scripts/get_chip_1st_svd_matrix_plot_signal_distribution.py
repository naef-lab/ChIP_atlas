import numpy as np
import pandas as pd
import os
import h5py
import matplotlib.pyplot as plt
import argparse
import multiprocessing

def parse_argument():
    parser = argparse.ArgumentParser(description='Make tensor of SVD 1st comp. of Chip signal')
    parser.add_argument('--infile_promoterome'
        ,required=True
        ,type=str
        ,help="Promoterome bed file")
    parser.add_argument('--infiles_svd'
        ,nargs='+'
        ,required=True
        ,type=str
        ,help="SVD files per TF")
    parser.add_argument('--outfile_tf_pos_prom'
        ,required=True
        ,type=str
        ,help="hdf5 with SVD 1st comp. tensor per TF, pos, prom")
    #parser.add_argument('--outfig'
    #    ,required=True
    #    ,type=str
    #    ,help="Figure with Chip signal distribution rel. to tss")
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['hg38','mm10']
        ,help="Genome version")
    parser.add_argument('--window_kb'
        ,required=True
        ,type=int
        ,help="Window size in kb")
    parser.add_argument('--threads'
        ,required=False
        ,type=int
        ,default=4
        ,help="Number of threads")
    return parser.parse_args()

def fill_in_tensor(infile):
    with h5py.File(infile,'r') as hf:
        U = hf['U'][:]
        S = hf['S'][:]
    return U*S


if __name__ == '__main__':

    # args
    args = parse_argument()
    
    promoterome = pd.read_csv(args.infile_promoterome,sep='\t')
    idx_minus_strand = (promoterome.strand=='-').values

    # get dims
    N_prom = promoterome.shape[0]
    with h5py.File(args.infiles_svd[0],'r') as hf:
        N_pos = hf['U'].shape[0]//N_prom
    N_tf = len(args.infiles_svd)
    
    # run loop in parralel
    #pool = multiprocessing.Pool(args.threads)
    #out = zip(*pool.map(fill_in_tensor, args.infiles_svd))
    #X = np.array(out).reshape([N_tf,N_pos,N_prom])
    #for t in range(N_tf):
    #    # flip promoters on minus strand
    #    X[t,:,idx_minus_strand] = X[t,:,idx_minus_strand][:,::-1]
    
    # initilize tf x pos x prom tensor
    X = np.zeros([N_tf,N_pos,N_prom])
    for t,infile in enumerate(args.infiles_svd):

        if t%20==0:
            print(np.round(t/N_tf,3))

        with h5py.File(infile,'r') as hf:
            u = hf['U'][:,0]
            s = hf['S'][0]
        tmp = np.reshape(u*s,[N_prom,N_pos]).T
        # flip promoters on minus strand
        tmp[:,idx_minus_strand] = tmp[:,idx_minus_strand][:,::-1]
        X[t,:,:] = tmp

    # save
    with h5py.File(args.outfile_tf_pos_prom,'w') as hf:
        hf.create_dataset('tf_pos_prom',data=X)

    #Chip_signal = X.mean(axis=2).mean(axis=0)
#
    ## get Chip signal cdf
    #cdf = np.cumsum(Chip_signal)/sum(Chip_signal)
#
    ## bin X accross position in equal intervals of the cdf
    #BIN_IDX = []
    #with h5py.File(args.outfile,'w') as hf:
    #    for N_bin in range(1,11):
    #        
    #        print(N_bin)
#
    #        if N_bin==1:
    #            X_binned = X.mean(axis=1)
    #        else:
    #            q = np.linspace(0,1,N_bin+1)[:-1]
    #            bin_idx  = np.argmin(np.abs( cdf[:,None] - q[None,:] ),axis=0)
#
    #            X_binned = np.zeros([N_tf,N_bin,N_prom])
    #            for i in range(N_bin-1):
    #                X_binned[:,i,:] = X[:,bin_idx[i]:bin_idx[i+1],:].mean(axis=1)
    #            X_binned[:,N_bin-1,:] = X[:,bin_idx[N_bin-1]:,:].mean(axis=1)
#
    #            X_binned = X_binned.reshape([N_tf*N_bin,N_prom])
    #            
    #            BIN_IDX.append(bin_idx)
#
    #        hf.create_dataset(str(N_bin),data=X_binned)

#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    
#    #x = np.arange(-int(window_kb)*1000,int(window_kb)*1000,20) + 10
#    x = np.linspace(-int(args.window_kb)*1000,int(args.window_kb)*1000,N_pos+1)
#    x = .5*(x[:-1]+x[1:])
#    bin_size = 2*int(args.window_kb)*1000/N_pos
#    ax.bar(x=x,height=Chip_signal,width=bin_size)
#    ms = np.linspace(8,4,9)
#    y = np.linspace(0,Chip_signal.min()/2,10)[1:]
#    for n in range(9):
#        ax.plot(x[BIN_IDX[n][1:]],np.repeat(y[n],n+1),'o',ms = ms[n],color=plt.cm.tab10(n+1))
#    ax.set_xlabel('Pos. rel. to tss')
#    ax.set_ylabel('Mean Chip signal per promoter and tf')
#    ax.grid('on')
#    fig.set_size_inches([8,6])
#    plt.tight_layout()
#    fig.savefig(args.outfig)


    