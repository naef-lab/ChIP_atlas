import pandas as pd
import numpy as np
import argparse
import pickle 

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName Synonyms dict.')
    parser.add_argument('--infile'
        ,required=True
        ,type=str)
    parser.add_argument('--outfile_table'
        ,required=True
        ,type=str)
    parser.add_argument('--outfile_dict'
        ,required=True
        ,type=str)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_argument()

    df = pd.read_csv(args.infile,sep='\t')
    df.replace({np.nan:'None'},inplace=True)

    # group by synonyms
    df_group = df.groupby(['Gene stable ID','Gene name']).agg({'Gene Synonym': lambda x : list(x)})
    df_group.loc[:,'Gene Synonym'] = df_group.loc[:,'Gene Synonym'].apply(lambda x: ','.join(x))
    df_group.to_csv(args.outfile_table,sep='\t')