import pyBigWig
import pandas as pd

infile='resources/experimentList_mm10_TFs.tab'
outfile='resources/to_dowload.txt'

chip_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)

with open(outfile,'w') as fout:
    for id in chip_tf.index:

        try:
            pyBigWig.open(f'resources/tracks/{id}.05.bb')
        except:
            print(id)
            fout.write(f'{id}.05.bb\n')
            
        try:
            pyBigWig.open(f'resources/tracks/{id}.bw')
        except:
            print(id)
            fout.write(f'{id}.bw\n')
            