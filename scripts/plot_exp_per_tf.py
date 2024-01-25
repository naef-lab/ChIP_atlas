import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    for genome in ['mm10','hg38']:
        df = pd.read_csv(f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab',sep='\t')

        # count exp per antigen
        N_exp_per_antigen = df.groupby('antigen').agg('count')['0'].values
        #N_exp_per_antigen = N_exp_per_antigen.sort_values('antigen')
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        bins = 2**( np.arange( np.ceil(np.log2(N_exp_per_antigen.max())) ) - .5 )
        ax.hist(N_exp_per_antigen,bins=bins)
        ax.set_xscale('log',base=2)
        ax.set_yscale('log')
        ax.set_xlabel(r'$log_{10}$ expermient per tf')
        ax.set_ylabel('nr. of tf')

        #N_exp_per_antigen.plot(kind='bar',ax=ax,width=1)
        fig.set_size_inches([6,3])
        plt.tight_layout()
        fig.savefig(f'results/fig/hist_experiment_per_tf_{genome}.pdf')
