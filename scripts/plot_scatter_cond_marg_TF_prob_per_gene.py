import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
     
     Genome = ['hg38','mm10']
     window = '2'
     N_bin = 1

     for genome in Genome:
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
        fig.savefig(f'results/fig/{genome}/scatter_Z_score_cond_marg_{N_bin}bin_{window}kb.pdf',dpi=300)
        plt.close()


