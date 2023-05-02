import pandas as pd

# input/output
infile='resources/experimentList_mm10_TFs.tab'
outfile = 'config/chip_seq.yaml'

# get tf:experiment dictionary
experiment_tf = pd.read_csv(infile,sep='\t',header=None,usecols=[0,3],index_col=0)
experiment_tf.columns = ['tf']
tf_dict = experiment_tf.groupby('tf').groups

# write configfile"
with open(outfile,'w') as fout:
    fout.write("promoter: 'resources/mm10_promoters_v2_pm1kb.gff'\n")
    fout.write("bin_size: 20\n")

    fout.write("TF_EXPERIMENT:\n")
    for tf in tf_dict.keys():
        fout.write(f"    {tf}:\n")
        for exp in tf_dict[tf]:
            fout.write(f"        - {exp}\n")