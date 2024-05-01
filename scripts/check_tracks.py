import pyBigWig
import pandas as pd
import argparse

#genome = 'hg38'
#genome = 'mm10'
#infile=f'resources/experimentList_{genome}_TFs_only_QC_filtered.tab'
#outfile=f'resources/to_dowload_{genome}.txt'

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName Synonyms dict.')
    parser.add_argument('--genome'
        ,required=True
        ,type=str)
    parser.add_argument('--infile'
        ,required=True
        ,type=str)
    parser.add_argument('--outfile'
        ,required=True
        ,type=str)
    
    return parser.parse_args()


if __name__ == '__main__':
    
    args = parse_argument()

    chip_tf = pd.read_csv(args.infile,sep='\t',usecols=[0,3],index_col=0)

    with open(args.outfile,'w') as fout:
        for id in chip_tf.index:

            try:
                pyBigWig.open(f'resources/tracks/{args.genome}/{id}.05.bb')
            except:
                print(id)
                fout.write(f'{id}.05.bb\n')
                
            try:
                pyBigWig.open(f'resources/tracks/{args.genome}/{id}.bw')
            except:
                print(id)
                fout.write(f'{id}.bw\n')
            