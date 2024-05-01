import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--infile'
        ,required=True
        ,type=str
        ,help="Experiments table")
    parser.add_argument('--genomes'
        ,required=True
        ,type=str
        ,nargs='+'
        ,help="Genomes to plot")
    parser.add_argument('--outfig_antigen_class_per_genome'
        ,required=True
        ,type=str
        ,help="Output figure pdf")
    parser.add_argument('--outfig_QC'
        ,required=True
        ,type=str
        ,help="Output figure pdf")
    parser.add_argument('--outfig_peaks_per_unique_mapped_reads'
        ,required=True
        ,type=str
        ,help="Output figure pdf")
    parser.add_argument('--th_reads'
        ,default=1000
        ,type=int)
    parser.add_argument('--th_mapped_reads'
        ,default=0.5
        ,type=float)
    parser.add_argument('--th_duplicates'
        ,default=0.5
        ,type=float)
    parser.add_argument('--th_peaks'
        ,default=10
        ,type=int)

    return parser.parse_args()


def plot_stats(chip,outfig):

    fig = plt.figure()
    cols = 1
    rows = 1
    ax = fig.add_subplot(rows,cols,1)

    pd.crosstab(chip.antigen_class,chip.genome).plot(kind='bar',ax=ax,width=0.7,cmap='tab20b')
    ax.set_yscale('log')
    ax.legend(loc='lower center',bbox_to_anchor=(.5, 1),title='genome',ncols=5)
    ax.set_ylabel('nr. of experiments')
        
    fig.set_size_inches([cols*8,rows*6])
    plt.tight_layout()
    fig.savefig(outfig)


if __name__ == '__main__':

    args = parse_argument()
    
    chip = pd.read_csv(args.infile,sep='\t',header=None,usecols=[0,1,2,3,4,5,6,7,8],index_col=0)
    chip.columns = ['genome','antigen_class','antigen','celltype_class','celltype','celltype_description','QC','title']
    # QC:
    #   ChIP/ATAC/DNase: nr. of reads, % mapped, % duplicates, # of peaks [Q < 1E-05]'
    #   Bisulfite-seq:   nr. of reads, % mapped, × coverage, # of hyper MR

    # plot stats
    outfig = 'results/fig/hist_antigen_class_per_genome.pdf'
    plot_stats(chip,args.outfig_antigen_class_per_genome)


    # get only mouse and human
    #chip = chip.loc[ (chip.genome=='hg19')|(chip.genome=='hg38')|(chip.genome=='mm10'),:]
    idx = np.any( np.array([(chip.genome==g).values for g in args.genomes]), axis=0 )
    chip = chip.loc[idx,:]

    # get only TFs and others
    chip = chip.loc[chip.antigen_class == 'TFs and others',:]

    QC = pd.DataFrame([[float(n) for n in qc.split(',')] for qc in chip.QC],columns=['n_reads','f_mapped','f_duplicates','n_peaks'],index=chip.index)
    QC.iloc[:,1] /= 100
    QC.iloc[:,2] /= 100
    QC['peaks_per_unique_mapped_reads'] = QC.n_peaks/(QC.n_reads * QC.f_mapped * (1-QC.f_duplicates))
    QC.loc[(QC.n_peaks==0) | (QC.f_mapped==0),'peaks_per_unique_mapped_reads'] = 0

    QC.iloc[:,0] = np.log10( QC.iloc[:,0] )
    QC.iloc[:,3] = np.log10( QC.iloc[:,3]+1 )
    QC.iloc[:,4] = np.log10( QC.iloc[:,4]+ 1e-6 )
    QC.columns = [r'$log_{10}$ nr. of reads', 'frac. mapped', 'frac. duplicates', r'$log_{10}$ nr. of peaks + 1',r'$log_{10}$ peaks/(uniq mapped read)']

    corr_QC = np.corrcoef( QC.values.T )

    chip = pd.concat([chip,QC],axis=1)

    scale = [['linear','log'],
             ['linear','log'],
             ['linear','linear'],
             ['linear','log'],
             ['linear','log']]
    
    thresholds = [np.log10(args.th_reads), args.th_mapped_reads, args.th_duplicates, np.log10(args.th_peaks)]

    Genome = chip.genome.unique()

    fig = plt.figure()
    cols = 2
    rows = 2

    for i,col in enumerate(QC.columns[:-1]):
        
        ax = fig.add_subplot(rows,cols,i+1)
        x = [ chip.loc[chip.genome == genome,col].values for genome in Genome ]
        if col==r'$log_{10}$ nr. of peaks + 1':
            for j in range(len(x)):
                x[j] = np.array(x[j])
                x[j][x[j]==0] = -1
                x[j] = list(x[j])

        ax.hist(x,bins=20,label=Genome)
        ylim = ax.get_ylim()
        th = [2*[thresholds[i]],ylim]
        ax.plot(th[0],th[1])
        ax.set_xscale(scale[i][0])
        ax.set_yscale(scale[i][1])
        ax.set_xlabel(QC.columns[i])
        ax.set_ylabel('nr. of exp.')
        if i==0:
            ax.legend(Genome)
    
    fig.set_size_inches([cols*4,rows*3])
    plt.tight_layout()
    fig.savefig(args.outfig_QC)

    i=4


    fig = plt.figure()
    cols = 1
    rows = 1
    ax = fig.add_subplot(rows,cols,1)
    
    col = QC.columns[-1]
    x = [ chip.loc[chip.genome == genome,col].values for genome in Genome ]

    ax.hist(x,bins=20,label=Genome)
    ax.set_xscale(scale[i][0])
    ax.set_yscale(scale[i][1])
    ax.set_xlabel(QC.columns[i])
    ax.set_ylabel('nr. of exp.')
    ax.legend(Genome)

    fig.set_size_inches([cols*6,rows*4])
    plt.tight_layout()
    fig.savefig(args.outfig_peaks_per_unique_mapped_reads)


    fig = plt.figure()
    cols = 2
    rows = 1
    i=0

    colors =  [ plt.cm.tab10(np.where(g==Genome)[0][0]) for g in chip.genome ]

    i+=1
    ax = fig.add_subplot(rows,cols,i)

    ax.scatter(x=chip['$log_{10}$ nr. of reads'],y=chip['$log_{10}$ nr. of peaks + 1'],c=colors,s=1,alpha=0.1,rasterized=True)
    ax.set_xlabel('$log_{10}$ nr. of reads')
    ax.set_ylabel('$log_{10}$ nr. of peaks + 1')
    ax.set_title(rf'$\rho = {np.round(corr_QC[0,3],2)}$')

    #col = r'$log_{10}$ peak_per_reads'
    #chip[col] = chip['$log_{10}$ nr. of peaks + 1'] - chip['$log_{10}$ nr. of reads']
    #x = [ chip.loc[chip.genome == genome,col].values for genome in Genome ]
    #ax.hist(x,bins=50,label=Genome)
    #QC.iloc[:,i].hist(bins=50,ax=ax)
    #ax.set_xscale('linear')
    #ax.set_yscale('log')
    #ax.set_xlabel(col)
    #ax.set_ylabel('nr. of exp.')
    #ax.legend(Genome)


    i+=1
    ax = fig.add_subplot(rows,cols,i)

    ax.scatter(x=chip['frac. mapped'],y=chip['$log_{10}$ nr. of peaks + 1'],c=colors,s=chip['$log_{10}$ nr. of reads'],alpha=0.1,rasterized=True)
    ax.set_ylabel('$log_{10}$ nr. of peaks + 1')
    ax.set_xlabel('frac. mapped')
    ax.set_title(rf'$\rho = {np.round(corr_QC[1,3],2)}$')
    
    fig.set_size_inches([cols*8,rows*6])
    plt.tight_layout()
    fig.savefig("results/fig/scatter_QC_vals",dpi=150)