import pandas as pd
import os
import pyBigWig

infile = "~/ChIPseq/resources/experimentList_mm10_TFs.tab"
Chip = pd.read_csv(infile,sep='\t',header=None,index_col=0,usecols=(0,3))
Chip.columns=['tf']

# sort by tf name
Chip = Chip.sort_values('tf',axis=0)
# take a small subset for testing
#Chip = Chip.loc[ (Chip.tf=='Clock') | (Chip.tf=='Arntl'), :]

outfile='mm10/genomes.txt'
outfile='hub.txt'
outfile='ChIP_Atlas.html'
outfile='url.txt'
outfile='mm10/trackDb.txt'
with open(outfile,'w', encoding="utf-8") as fout:
    for tf in Chip.tf.unique():
        print(tf)

        # BigWig tracks

        fout.write(f"track ChIP_{tf}\n")
        fout.write("type bigWig\n")
        fout.write("container multiWig\n")
        fout.write(f"shortLabel {tf}\n")
        fout.write(f"longLabel ChIP-seq {tf}\n")
        fout.write("visibility hide\n")
        fout.write("aggregate transparentOverlay\n")
        fout.write("showSubtrackColorOnUi on\n")
        fout.write("maxHeightPixels 500:100:8\n")
        fout.write("autoScale on\n")
        fout.write("\n")

        i = 0

        n = sum(Chip.tf==tf)
        G = min(n+159,240)
        for idx in Chip.loc[Chip.tf==tf].index:
            if os.path.exists(f'tracks/{idx}.bw'):
                try:
                    pyBigWig.open(f'tracks/{idx}.bw')
                except:
                    print('file corrupted')
                finally:
                    i+=1

                    fout.write(f"\ttrack ChIP_{tf}_{i}\n")
                    fout.write(f"\tbigDataUrl http://upnaesrv1.epfl.ch/ChIP_Atlas/tracks/{idx}.bw\n")
                    fout.write(f"\tshortLabel {tf}_{i}\n")
                    fout.write(f"\tlongLabel Chip_{tf}_{i}\n")
                    fout.write(f"\tparent ChIP_{tf}\n")
                    fout.write("\ttype bigWig\n")
                    fout.write(f"\tcolor 255,{G},255\n")
                    fout.write("\n")

        fout.write("###################\n\n")
