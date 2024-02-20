import numpy as np
import pandas as pd
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make bed files for TF peaks')
    parser.add_argument('--genome', type=str, default='mm10', help='Genome')
    return parser.parse_args()

if __name__ == '__main__':

    CHR = [f'chr{c}' for c in list(range(1,23))+['X','Y','MT']]

    genome = parse_args().genome

    experimentList = pd.read_csv(f'/bigdata/jbreda/ChIP/resources/experimentList_{genome}_TFs_only_QC_filtered.tab',sep='\t',index_col=0)

    in_folder = f'/bigdata/jbreda/ChIP/resources/tracks/{genome}'
    out_folder = f'/bigdata/jbreda/ChIP/results/{genome}/TF_peaks'

    for tf in experimentList.antigen.unique():
        print(tf)

        bed_file = f'{out_folder}/{tf}.bed'
        if os.path.exists(bed_file):
            os.system(f'rm {bed_file}')
        os.system(f'touch {bed_file}')
        for exp_id in experimentList.loc[experimentList.antigen == tf,:].index:
            os.system(f'bigBedToBed {in_folder}/{exp_id}.05.bb {in_folder}/tmp.bed; cat {in_folder}/tmp.bed >> {bed_file}')
            os.system(f'rm {in_folder}/tmp.bed')
        
        bed = pd.read_csv(bed_file,sep='\t',header=None)
        bed.columns = ['chrom','chromStart','chromEnd','score']

        # rename chrM to chrMT
        bed['chrom'] = bed['chrom'].replace('chrM','chrMT')

        # keep only chromosomes in CHR
        bed = bed[ bed.chrom.isin(CHR) ]


        #score = bed['score'].values
        #score = (score - score.min()) / (score.max() - score.min()) 
        #rgb = np.array([1-score,score,0]).T
        rgb = np.zeros([bed.shape[0],3])
    
        # Add columns
        bed['name'] = tf
        bed['strand'] = '.'
        bed['thickStart'] = bed['chromStart']
        bed['thickEnd'] = bed['chromEnd']
        bed['itemRgb'] = [','.join(c) for c in (255*rgb).astype(int).astype(str)]
        bed['blockCount'] = 1
        bed['blockSizes'] = bed['chromEnd'] - bed['chromStart']
        bed['blockStarts'] = 0

        # clip score to 1000
        bed['score'] = np.clip(bed['score'],0,1000)


        # reorder columns
        bed_cols = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
        bed = bed[bed_cols]
        
        # save bed and bigbed
        chrom_sizes = {'mm10': '/bigdata/jbreda/genome/mm10/Mus_musculus.GRCm38.dna.primary_assembly.chrom.sizes',
                       'hg38': '/bigdata/jbreda/genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes'}
        out_bed = f'{out_folder}/{tf}.bed'
        out_bb = f'{out_folder}/{tf}.bb'
        bed.to_csv(out_bed,sep='\t',index=False,header=False)
        os.system(f'bedSort {out_bed} {out_bed}')
        os.system(f'bedToBigBed {out_bed} {chrom_sizes[genome]} {out_bb}')