import numpy as np
import pandas as pd
import argparse
import pyBigWig
import h5py

def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP Peaks for tf in tensor (N_promoters x N_positions X N_experiments).')
    parser.add_argument('--promoterome'
        ,required=True
        ,type=str
        ,help="Promoterome file")
    parser.add_argument('--tf'
        ,required=True
        ,type=str
        ,help="TF to get experiments")
    parser.add_argument('--genome'
        ,required=True
        ,type=str
        ,choices=['mm10','hg38']
        ,help="genome version")
    parser.add_argument('--infiles_tf'
        ,required=True
        ,type=str
        ,nargs='+'
        ,help="experiments peak files")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output table with peaks")

    return parser.parse_args()

if __name__ == '__main__':

    # Make a boolean tensor of shape (N_prom,N_pos,N_exp)
    # where N_prom is the number of promoters, N_pos is the number of base pairs in the window, and N_exp is the number of experiments
    # The tensor is True if a peak is found in the promoterome

    args = parse_argument()

    # load promoters
    promoterome = pd.read_csv(args.promoterome ,sep='\t')
    CHR = promoterome.chr.unique()

    # make peaks table from bb files
    peaks_table = pd.DataFrame(columns=['chr','start','end','gene','score','strand','promoter_id','antigen','exp_id','rel_start','rel_end'])
    dtypes = {'chr':'str',
                'start':'int',
                'end':'int',
                'gene':'str',
                'score':'float',
                'strand':'str',
                'promoter_id':'str',
                'antigen':'str',
                'exp_id':'str',
                'rel_start':'int',
                'rel_end':'int'}
    peaks_table = peaks_table.astype(dtypes)

    for infile in args.infiles_tf:
        try:
            bb = pyBigWig.open(infile)
        except:
            print(f"Cannot open {infile}")
            continue

        exp_id = infile.split('/')[-1].split('.')[0]
        
        for chr in CHR:
            if chr in bb.chroms():
                # get promoterome coordinates
                my_prom = promoterome[promoterome.chr==chr]
                x = my_prom.loc[:,['start','end']].values

                # get peaks coordinates
                my_peaks = pd.DataFrame(bb.entries(chr,0,bb.chroms(chr)),columns=['start','end','score'])
                y = my_peaks.loc[:,['start','end']].values

                # get overlap between peaks and promoterome
                Prom_Peak_overlap = (x[:,0][:,None] <= y[:,1][None,:]) & (x[:,1][:,None] >= y[:,0][None,:])
                [idx_prom,idx_peak] = np.where(Prom_Peak_overlap)

                if len(idx_prom) == 0:
                    continue

                # get position of peaks relative to tss
                tss = (my_prom.iloc[idx_prom,:].start.values + my_prom.iloc[idx_prom,:].end.values)//2
                rel_start = my_peaks.iloc[idx_peak,0].values - tss
                rel_end = my_peaks.iloc[idx_peak,1].values - tss
                # Now fill in with overlaping peaks
                my_dict = {'chr':chr,
                           'start':my_peaks.iloc[idx_peak,0].values,
                           'end':my_peaks.iloc[idx_peak,1].values,
                           'gene':my_prom.iloc[idx_prom,:].gene.values,
                           'score':my_peaks.iloc[idx_peak,2].values,
                           'strand':my_prom.iloc[idx_prom,:].strand.values,
                           'promoter_id': my_prom.iloc[idx_prom,:].id.values,
                           'antigen':args.tf,
                           'exp_id':exp_id,
                           'rel_start':rel_start,
                           'rel_end':rel_end}

                # append to peaks table
                peaks_table = pd.concat([peaks_table,pd.DataFrame(my_dict).astype(dtypes)],axis=0)

    # save
    peaks_table.to_csv(args.outfile,sep='\t',index=False)


