import pandas as pd
import argparse
import numpy as np
import pickle

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName Synonyms dict.')
    parser.add_argument('--infile_chip'
        ,required=True
        ,type=str)
    parser.add_argument('--infile_tfs'
        ,required=True
        ,type=str)
    parser.add_argument('--geneid_genename_synonym_table'
        ,required=True
        ,type=str)
    parser.add_argument('--synonym_genename_dict'
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
    parser.add_argument('--th_exp_per_tf'
        ,default=250
        ,type=int)
    
    return parser.parse_args()

def get_experiment_table(args):
    
    # get full experiment table
    chip = pd.read_csv(args.infile_chip,sep='\t',header=None,usecols=[0,1,2,3,4,5,6,7,8],index_col=0)
    chip.columns = ['genome','antigen_class','antigen','celltype_class','celltype','celltype_description','QC','title']
    
    # get only TFs and others from genome
    chip = chip[ (chip.genome == args.genome) & (chip.antigen_class=='TFs and others') ]

    # parse QC column add to chip table
    QC = pd.DataFrame([[float(n) for n in qc.split(',')] for qc in chip.QC],columns=['n_reads','f_mapped','f_duplicates','n_peaks'],index=chip.index)
    QC.iloc[:,1] /= 100
    QC.iloc[:,2] /= 100
    chip = pd.concat([chip,QC],axis=1)

    # get peaks per unique mapped reads
    chip.loc[:,'n_peaks_per_unique_mapped_reads'] = chip.n_peaks/(chip.f_mapped*chip.n_reads*(1-chip.f_duplicates))
    n_tot = chip.shape[0]
    print(f'{args.genome} - {n_tot} experiments')

    # apply thresholds
    idx_out = list( chip[(chip['n_reads']     < args.th_reads) | 
                         (chip['f_mapped']    < args.th_mapped_reads) |
                         (chip['f_duplicates']> args.th_duplicates) |
                         (chip['n_peaks']     < args.th_peaks) ].index )
    
    # For TFs with more than th experiments, keep the th with the highest n_peaks_per_unique_mapped_reads
    exp_per_tf = chip.groupby('antigen')['antigen'].aggregate('count')
    for tf in exp_per_tf.loc[exp_per_tf > args.th_exp_per_tf].index:
        idx_out.extend( list( chip.loc[chip.antigen==tf].sort_values('n_peaks_per_unique_mapped_reads',ascending=False).index[args.th_exp_per_tf:] ) )
    chip.drop(idx_out,inplace=True)

    print(f'{args.genome} - {chip.shape[0]/n_tot} passed QC')
    return chip

if __name__ == '__main__':
    
    args = parse_argument()

    # get experiment table with QC and filtered fur up to th_exp_per_tf experiments per TF
    chip = get_experiment_table(args)

    # load GeneID GeneName Synonym table gene name list
    Gene_id_name_syn = pd.read_csv(args.geneid_genename_synonym_table,sep='\t')

    # get all gene names (Gene symbol)
    GeneName = set(Gene_id_name_syn['Gene name'])

    # remove genes with no synonyms, duplicates and with same name as synonym
    Gene_id_name_syn.replace({np.nan:'None'},inplace=True)
    Gene_id_name_syn.drop(Gene_id_name_syn.index[Gene_id_name_syn['Gene Synonym']=="None"],inplace=True)
    Gene_id_name_syn.drop_duplicates(inplace=True)
    Gene_id_name_syn.drop( Gene_id_name_syn.loc[Gene_id_name_syn['Gene Synonym'] == Gene_id_name_syn['Gene name']].index, inplace=True)

    # load synonym to gene name dictionary
    with open(args.synonym_genename_dict, 'rb') as f:
        Synonym_2_GeneName = pickle.load(f)

    # get antigen gene names to change
    to_rename = []
    not_found = []
    for g in chip.antigen.unique():
        if g in GeneName:
            continue
        else:
            if g in Synonym_2_GeneName.keys():
                to_rename.append([g,Synonym_2_GeneName[g]])
            else:
                not_found.append(g)
    to_rename = dict(zip(np.array(to_rename)[:,0],np.array(to_rename)[:,1]))

    # print to rename and not found
    print(f'{args.genome} - {len(to_rename)} antigens to rename')
    print(f'{args.genome} - {len(not_found)} antigens not found')

    # change gene names
    for g in to_rename.keys():
        idx = chip[ chip.antigen==g ].index
        for i in idx:
            chip.at[i,'antigen'] = to_rename[chip.at[i,'antigen']]

    # load TF list
    with open(args.infile_tfs,'r') as f:
        TFs = np.array( [tf.strip() for tf in f.readlines()] )
    
    # get tf names to change
    to_rename = []
    not_found = []
    for g in TFs:
        if g in GeneName:
            continue
        else:
            if g in Synonym_2_GeneName.keys():
                to_rename.append([g,Synonym_2_GeneName[g]])
            else:
                not_found.append(g)
    to_rename = dict(zip(np.array(to_rename)[:,0],np.array(to_rename)[:,1]))

    # change tf names
    for g in to_rename.keys():
        i = np.where( TFs == g )[0]
        TFs[i] = to_rename[g]

    # print to rename and not found
    print(f'{args.genome} - {len(to_rename)} TFs to rename')
    print(f'{args.genome} - {len(not_found)} TFs not found')

    # Keep only antigens that are in TF list
    tf_id = []
    tf_out = []
    for id in chip.index:
        antigen = chip.at[id,'antigen']
        # if antigen is in TF list add id
        if antigen in TFs:
            tf_id.append(id)
        # if antigen is a synonym of 
        elif antigen in Gene_id_name_syn['Gene Synonym'].values:
            gene_name = Gene_id_name_syn.loc[ Gene_id_name_syn['Gene Synonym']==antigen, 'Gene name'].values
            if any([g in TFs for g in gene_name]):
                tf_id.append(id)
        else:
            tf_out.append(antigen)

    # print kept ratio
    print(f'{args.genome} - {len(tf_id)/chip.shape[0]} in TF list')
    chip = chip.loc[tf_id,:]
    print(f'{args.genome} - {len(chip.antigen.unique())} unique TF')
    
    # write output table
    chip.to_csv(args.outfile,sep='\t')