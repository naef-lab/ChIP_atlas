import pandas as pd
import numpy as np
import pickle
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName TF list.')
    parser.add_argument('--geneid_genename_synonym_table'
        ,required=True
        ,type=str)
    parser.add_argument('--synonym_genename_dict'
        ,required=True
        ,type=str)
    parser.add_argument('--infile_tf_list'
        ,required=True
        ,type=str)
    parser.add_argument('--infiles_GO_terms'
        ,nargs='+'
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

    # ENSID to gene name
    df = pd.read_csv(args.geneid_genename_synonym_table,sep='\t',usecols=[0,1])
    df.replace({np.nan:'None'},inplace=True)
    ID_2_GeneName = dict(zip(df['Gene stable ID'],df['Gene name']))

    # load synonym to gene name dictionary
    with open(args.synonym_genename_dict, 'rb') as f:
        Synonym_2_GeneName = pickle.load(f)

    # load GO terms
    go_genes = set()
    for f in args.infiles_GO_terms:
        with open(f,'r') as fin:
            for line in fin:
                go_genes.add(line.strip())

    # correct gene names
    to_add = set()
    to_remove = set()
    for g in go_genes:
        if g in ID_2_GeneName.values():
            continue
        else:
            if g in Synonym_2_GeneName.keys():
                to_add.add(Synonym_2_GeneName[g])
            else:
                to_remove.add(g)

    go_genes = go_genes.difference(to_remove).union(to_add)

    # get curated TF list
    # mouse
    if args.genome == 'mm10':

        TFs = set()
        with open(args.infile_tf_list,'r') as fin:
            for line in fin:
                tf = line.strip()
                if tf in ID_2_GeneName:
                    TFs.add(ID_2_GeneName[tf])

    # human
    if args.genome == 'hg38':

        df = pd.read_excel(args.infile_tf_list, sheet_name="Table S1. Related to Figure 1B", skiprows=1,index_col=0)
        # rename column 3:
        rename = {"Unnamed: 3": "Is TF"}
        df = df.rename(columns=rename)

        TFs = set( np.sort( df.loc[df.loc[:,"Is TF"] == "Yes",:].Name.unique() ) )

        # correct gene names
        to_add = set()
        to_remove = set()
        for g in TFs:
            if g in ID_2_GeneName.values():
                continue
            else:
                if g in Synonym_2_GeneName.keys():
                    to_add.add(Synonym_2_GeneName[g])
                else:
                    to_remove.add(g)
        TFs = TFs.difference(to_remove).union(to_add)

    # join TFs and go_genes
    TFs = go_genes.union(set(TFs))
    
    with open(args.outfile,'w') as fout:
        for tf in TFs:
            fout.write(f'{tf}\n')