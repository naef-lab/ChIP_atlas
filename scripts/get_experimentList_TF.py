import pandas as pd
import argparse
import numpy as np

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName Synonyms dict.')
    parser.add_argument('--infile_chip'
        ,required=True
        ,type=str)
    parser.add_argument('--infile_tfs'
        ,required=True
        ,type=str)
    parser.add_argument('--infile_gene_dict'
        ,required=True
        ,type=str)
    parser.add_argument('--genome'
        ,required=True
        ,type=str)
    parser.add_argument('--outfile'
        ,required=True
        ,type=str)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_argument()

    # get experiment table
    chip = pd.read_csv(args.infile_chip,sep='\t',header=None,usecols=[0,1,2,3],index_col=0)
    chip.columns = ['genome','Antigen_class','Antigen']
    chip = chip[ (chip.genome == args.genome) & (chip.Antigen_class=='TFs and others') ]

    # get TF list
    TFs = pd.read_csv(args.infile_tfs,sep='\t',header=None,usecols=[0,1])
    TFs.columns = ['GeneID','GeneName']

    # get Gene ID/name/synonyms table
    Gene_id_name_syn = pd.read_csv(args.infile_gene_dict,sep='\t')
    
    # keep only Gene id that are in TFs list
    idx_in = [Gene_id_name_syn.at[i,'Gene stable ID'] in TFs.GeneID.values for i in Gene_id_name_syn.index]
    Gene_id_name_syn = Gene_id_name_syn.loc[idx_in,:]
    
    # change nans -> 'None', remove None
    Gene_id_name_syn.replace({np.nan:'None'},inplace=True)
    idx_out = Gene_id_name_syn['Gene Synonym']=='None'
    Gene_id_name_syn = Gene_id_name_syn.loc[~idx_out,:]
    
    # remove and synonym same as name
    idx_out = Gene_id_name_syn['Gene name'] == Gene_id_name_syn['Gene Synonym']
    Gene_id_name_syn = Gene_id_name_syn.loc[~idx_out,:]

    
    Gene_id_name_syn.drop_duplicates(inplace=True)


    tf_id = []
    for id in chip.index:
        antigen = chip.at[id,'Antigen']
        # if antigen is in TF list add id
        if antigen in TFs.GeneName.values:
            tf_id.append(id)
        # if antigen is a synonym of 
        elif antigen in Gene_id_name_syn['Gene Synonym'].values:
            gene_name = Gene_id_name_syn.loc[ Gene_id_name_syn['Gene Synonym']==antigen, 'Gene name'].values
            if gene_name in TFs.GeneName.values:
                tf_id.append(id)

    print(f'{len(tf_id)} / {chip.shape[0]} genes in TF list')
    chip = chip.loc[tf_id,:]

    print(f'{len(chip.Antigen.unique())} unique TF')
    
    chip.to_csv(args.outfile,sep='\t')