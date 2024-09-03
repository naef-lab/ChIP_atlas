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

if __name__ == '__main__':
    
    args = parse_argument()

    # get experiment table
    chip = pd.read_csv(args.infile_chip,sep='\t',header=None,usecols=[0,1,2,3,4,5,6,7,8],index_col=0)
    chip.columns = ['genome','antigen_class','antigen','celltype_class','celltype','celltype_description','QC','title']
    
    # get only TFs and others from genome
    chip = chip[ (chip.genome == args.genome) & (chip.antigen_class=='TFs and others') ]

    # load GeneID GeneName Synonym table gene name list
    Gene_id_name_syn = pd.read_csv(args.geneid_genename_synonym_table,sep='\t')
    Gene_id_name_syn.replace({np.nan:'None'},inplace=True)
    Gene_id_name_syn.drop(Gene_id_name_syn.index[Gene_id_name_syn['Gene Synonym']=="None"],inplace=True)
    Gene_id_name_syn.drop_duplicates(inplace=True)
    Gene_id_name_syn.drop( Gene_id_name_syn.loc[Gene_id_name_syn['Gene Synonym'] == Gene_id_name_syn['Gene name']].index, inplace=True)

    GeneName = set(Gene_id_name_syn['Gene name'])

    # load synonym to gene name dictionary
    with open(args.synonym_genename_dict, 'rb') as f:
        Synonym_2_GeneName = pickle.load(f)

    # parse QC column add to chip table
    QC = pd.DataFrame([[float(n) for n in qc.split(',')] for qc in chip.QC],columns=['n_reads','f_mapped','f_duplicates','n_peaks'],index=chip.index)
    QC.iloc[:,1] /= 100
    QC.iloc[:,2] /= 100
    chip = pd.concat([chip,QC],axis=1)

    # get peaks per unique mapped reads
    chip.loc[:,'n_peaks_per_unique_mapped_reads'] = chip.n_peaks/(chip.f_mapped*chip.n_reads*(1-chip.f_duplicates))
    n_tot = chip.shape[0]
    print(f'{args.genome}: {n_tot} experiments')

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

    print(f'{args.genome}: {chip.shape[0]/n_tot} passed QC')

    # get TF list
    with open(args.infile_tfs,'r') as f:
        TFs = [tf.strip() for tf in f.readlines()]

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

    # change gene names
    for g in to_rename.keys():
        idx = chip[ chip.antigen==g ].index
        for i in idx:
            chip.at[i,'antigen'] = to_rename[chip.at[i,'antigen']]

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
    print(f'{args.genome}: {len(tf_id)/chip.shape[0]} in TF list')
    chip = chip.loc[tf_id,:]
    print(f'{args.genome}: {len(chip.antigen.unique())} unique TF')

    # Take top experiments according to n_peaks_per_unique_mapped_reads for TFs with more than 250 experiments
    n_tot = chip.shape[0]
    TF,N = np.unique(chip.antigen,return_counts=True)
    for my_tf in TF[N>args.th_exp_per_tf]:
        print(my_tf)
        my_chip = chip.loc[chip.antigen==my_tf,['n_peaks_per_unique_mapped_reads']]
        th = np.sort(my_chip.n_peaks_per_unique_mapped_reads.values)[-args.th_exp_per_tf]

        idx_out = my_chip.loc[my_chip.n_peaks_per_unique_mapped_reads<th,:].index
        chip.drop(index=idx_out,inplace=True)
    
    print(f'{args.genome}: {chip.shape[0]/n_tot} atfer shortening max nr. exp per TF')
    
    # write output table
    chip.to_csv(args.outfile,sep='\t')