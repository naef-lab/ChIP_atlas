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
    df.drop(df.index[df['Gene Synonym']=="None"],inplace=True)

    # group by synonyms
    #df_group = df.groupby(['Gene stable ID','Gene name']).agg({'Gene Synonym': lambda x : list(x)})
    df_group = df.groupby(['Gene name']).agg({'Gene Synonym': lambda x : set(x)})

    # save table in text format
    df_group.loc[:,'Gene Synonym'].apply(lambda x: ','.join(x)).to_csv(args.outfile_table,sep='\t')

    # remove synonyms with more than one gene name
    all_synonyms = []
    for syn in df_group['Gene Synonym']:
        all_synonyms.extend(list(syn))
    u, c = np.unique(np.array(all_synonyms), return_counts=True)
    to_remove = set(u[c>1])
    for g in df_group.index:
        df_group.at[g,'Gene Synonym'] -= to_remove
    
    # make dict
    my_dict = {}
    for g in df_group.index:
        for syn in df_group.at[g,'Gene Synonym']:
            my_dict[syn] = g

    # save dict
    with open(args.outfile_dict, 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)