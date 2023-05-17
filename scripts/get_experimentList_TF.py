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

if __name__ == '__main__':
    
    args = parse_argument()

    # get experiment table
    chip = pd.read_csv(args.infile_chip,sep='\t',header=None,usecols=[0,1,2,3,4,5,6,7,8],index_col=0)
    chip.columns = ['genome','antigen_class','antigen','celltype_class','celltype','celltype_description','QC','title']
    
    # get only TFs and others from genome
    chip = chip[ (chip.genome == args.genome) & (chip.antigen_class=='TFs and others') ]
    
    # parse QC column add to chip table
    QC = pd.DataFrame([[float(n) for n in qc.split(',')] for qc in chip.QC],columns=['n_reads','f_mapped','f_duplicates','n_peaks'],index=chip.index)
    QC.iloc[:,1] /= 100
    QC.iloc[:,2] /= 100
    chip = pd.concat([chip,QC],axis=1)

    n_tot = chip.shape[0]
    # apply thresholds
    chip = chip[(chip['n_reads']     > args.th_reads) & 
                (chip['f_mapped']    > args.th_mapped_reads) &
                (chip['f_duplicates']< args.th_duplicates) &
                (chip['n_peaks']     > args.th_peaks) ]
    
    print(f'{args.genome}: {chip.shape[0]/n_tot} passed QC')

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

    # remove duplicates
    Gene_id_name_syn.drop_duplicates(inplace=True)

    # Keep only antigens that are in TF list or which have synonym in TF list
    tf_id = []
    for id in chip.index:
        antigen = chip.at[id,'antigen']
        # if antigen is in TF list add id
        if antigen in TFs.GeneName.values:
            tf_id.append(id)
        # if antigen is a synonym of 
        elif antigen in Gene_id_name_syn['Gene Synonym'].values:
            gene_name = Gene_id_name_syn.loc[ Gene_id_name_syn['Gene Synonym']==antigen, 'Gene name'].values
            if gene_name in TFs.GeneName.values:
                tf_id.append(id)

    # print kept ratio
    print(f'{args.genome}: {len(tf_id)/chip.shape[0]} in TF list')
    chip = chip.loc[tf_id,:]
    print(f'{args.genome}: {len(chip.antigen.unique())} unique TF')
    
    chip.to_csv(args.outfile,sep='\t')