import pyBigWig
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make bed files for TF peaks')
    parser.add_argument('--genome', type=str, default='mm10', help='Genome')
    parser.add_argument('--promoterome', type=str, help='Promoterome')
    parser.add_argument('--out_file', type=str, help='Output file')
    return parser.parse_args()

if __name__ == '__main__':

    args = parse_args()

    # get experiment table and promoterome
    experiments = pd.read_csv(f'resources/experimentList_{args.genome}_TFs_only_QC_filtered.tab',sep = '\t',index_col = 0)
    promoterome = pd.read_csv(args.promoterome,sep='\t')
    if args.genome == 'mm10':
        CHR = [f'chr{i+1}' for i in range(19)] + ['chrX','chrY','chrM']
    elif args.genome == 'hg38':
        CHR = [f'chr{i+1}' for i in range(22)] + ['chrX','chrY','chrM']

    # get peaks
    peaks_table = pd.DataFrame(columns=['exp_id','chr','start','end','score'])
    peak_infiles = [f'resources/tracks/{args.genome}/{id}.05.bb' for id in experiments.index]
    for infile in peak_infiles:
        print(infile)
        bb = pyBigWig.open(infile)
        id = infile.split('/')[-1].split('.')[0]
        for c in CHR:
            if c in bb.chroms():
                pks = pd.DataFrame(bb.entries(c,0,bb.chroms(c)),columns=['start','end','score'])
                pks['chr'] = c
                pks['exp_id'] = id
                peaks_table = pd.concat([peaks_table,pks],axis=0)
    peaks_table = peaks_table.reset_index(drop=True)

    # write peaks to bed file
    peaks_table.to_csv(args.out_file,sep='\t',index=False,header=True)

    