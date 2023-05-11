import pandas as pd
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Make GeneID GeneName TF list.')
    parser.add_argument('--infile_tf'
        ,required=True
        ,type=str)
    parser.add_argument('--infile_gene_dict'
        ,required=True
        ,type=str)
    parser.add_argument('--outfile'
        ,required=True
        ,type=str)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_argument()

    Gene_dict = pd.read_csv(args.infile_gene_dict,sep='\t',index_col=0)
    geneid2name = Gene_dict['Gene name'].to_dict()

    TFs = pd.read_csv(args.infile_tf,sep='\t',header=None,index_col=0)
    TFs['GeneName'] = ""

    for id in TFs.index:
        if id in geneid2name:
            TFs.at[id,'GeneName'] = geneid2name[id]
        else:
            TFs.at[id,'GeneName'] = 'unknown'
            print(id,' not found')
    
    TFs.to_csv(args.outfile,sep='\t',header=False, index=True)