import pandas as pd
import os

if __name__ == '__main__':

    # config file
    outfile = 'config/chip_seq.yaml'
    template = 'config/template_chip_seq.yaml'

    GENOME = ['mm10','hg19','hg38']
    tab = "  "

    # copy template
    os.system(f'cp {template} {outfile}')

    # write configfile
    with open(outfile,'a') as fout:

        fout.write("TF_EXPERIMENT:\n")

        for genome in GENOME:
            print(genome)

            fout.write(f"{tab}{genome}:\n")

            # input/output
            infile=f'resources/experimentList_{genome}_TFs_only.tab'

            # get tf:experiment dictionary
            experiment_tf = pd.read_csv(infile,sep='\t',usecols=[0,3],index_col=0)
            experiment_tf.columns = ['tf']
            tf_dict = experiment_tf.groupby('tf').groups
            
            for tf in tf_dict.keys():
                fout.write(f"{tab}{tab}{tf}:\n")
                for exp in tf_dict[tf]:
                    fout.write(f"{tab}{tab}{tab}- {exp}\n")