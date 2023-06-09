import pyBigWig
import pandas as pd

genome = 'hg38'

infile=f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab'
outfile=f'resources/to_dowload_{genome}.txt'

chip_tf = pd.read_csv(infile,sep='\t',usecols=[0,3],index_col=0)

with open(outfile,'w') as fout:
    for id in chip_tf.index:

        try:
            pyBigWig.open(f'resources/tracks/{genome}/{id}.05.bb')
        except:
            print(id)
            fout.write(f'{id}.05.bb\n')
            
        try:
            pyBigWig.open(f'resources/tracks/{genome}/{id}.bw')
        except:
            print(id)
            fout.write(f'{id}.bw\n')
            