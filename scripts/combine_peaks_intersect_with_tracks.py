import numpy as np
import pandas as pd
import pyBigWig
import os

CHR=[f'chr{c}' for c in list(range(1,20))+['X','Y']]

if __name__ == "__main__":
    
    # load chip experiments table
    chip_tf = pd.read_csv('resources/experimentList_mm10_TFs.tab',sep='\t',header=None,usecols=[0,3],index_col=0)
    chip_tf.columns = ['tf']

    # get unique tf
    TF = chip_tf.tf.unique()

    for tf in TF:
        print(tf)
        
        # get indices of tf
        idx = chip_tf.loc[chip_tf.tf == tf,:].index

        # combine if > 1 experiment
        if len(idx) > 1:
            bb_out = pyBigWig.open(f'resources/merged_peaks/{tf}.05.bb','w')
            for i in idx:
                bb = pyBigWig.open(f'resources/tracks/{i}.05.bb')
        else:
            os.system(f'cp resources/tracks/{i}.05.bb resources/merged_peaks/{tf}.05.bb')
