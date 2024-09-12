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
    parser.add_argument('--window_kb'
        ,default=2
        ,type=int
        ,choices=[1,2,3,5]
        ,help="window size (in kb)")
    parser.add_argument('--infiles_tf'
        ,required=True
        ,type=str
        ,nargs='+'
        ,help="experiments peak files")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output hfd5 file with tensor")

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
    peaks_table = pd.DataFrame(columns=['exp_id','chr','start','end','score'])
    bb_files = []
    for infile in args.infiles_tf:
        try:
            bb = pyBigWig.open(infile)
        except:
            print(f"Cannot open {infile}")
            continue

        bb_files.append(infile)

        id = infile.split('/')[-1].split('.')[0]
        for chr in CHR:
            if chr in bb.chroms():
                p = pd.DataFrame(bb.entries(chr,0,bb.chroms(chr)),columns=['start','end','score'])
                p['chr'] = chr
                p['exp_id'] = id
                peaks_table = pd.concat([peaks_table,p],axis=0)
    peaks_table = peaks_table.reset_index(drop=True)

    # get tensor dimension and initialize
    N_exp = len(bb_files)
    N_prom = promoterome.shape[0]
    N_pos = promoterome.at[0,'end'] - promoterome.at[0,'start']
    # tensor of shape (N_prom,N_pos,N_exp) X[i,j,k] is True if peak k is found in promoter i at position j
    TF_peaks = np.zeros([N_prom,N_pos,N_exp],dtype=bool)
    
    exp_ids = np.array( [exp.split('/')[-1].split('.')[0] for exp in bb_files] )
    x = promoterome.loc[:,['chr','start','end']].values

    for k, exp_id in enumerate(exp_ids):
        
        # get peaks for experiment
        peak_idx = peaks_table[peaks_table.exp_id==exp_id].index
        y = peaks_table.loc[peak_idx,['chr','start','end']].values

        # get overlap between peaks and promoterome
        Prom_Peak_overlap = (x[:,0][:,None] == y[:,0][None,:]) & (x[:,1][:,None] <= y[:,2][None,:]) & (x[:,2][:,None] >= y[:,1][None,:])
        [idx_prom,idx_peak] = np.where(Prom_Peak_overlap)

        # Now fill in the tensor with the overlaping peaks
        for i,p in zip(idx_prom,idx_peak):
            o_start = peaks_table.loc[peak_idx[p],'start'] - promoterome.loc[i,'start']
            o_end = peaks_table.loc[peak_idx[p],'end'] - promoterome.loc[i,'start']
            j = range(max(o_start,0),min(o_end,N_pos))
            TF_peaks[i,j,k] = True

    # write tensor in out h5py file after removing failed expeiment
    with h5py.File(args.outfile,'w') as hf:
        d = hf.create_dataset('peak_prom_pos_exp',data=TF_peaks)
        # add metadata
        d.attrs['experiment_id'] = ','.join(exp_ids)
        d.attrs['TF'] = args.tf
        d.attrs['genome'] = args.genome
        d.attrs['window_kb'] = args.window_kb
        d.attrs['promoterome'] = args.promoterome
