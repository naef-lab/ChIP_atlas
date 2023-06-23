import pandas as pd
import os
import pyBigWig


if __name__ == '__main__':

    Genomes = ['mm10','hg38']
    outfolder = 'track_hubs'

    # make track hub
    outfile=f'{outfolder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("hub ChIP_Atlas\n")
        fout.write("shortLabel ChIP_atlas\n")
        fout.write("longLabel A ChIP-seq atlas of TFs in mouse and human tissues\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write("descriptionUrl ChIP_Altas.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{outfolder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        for genome in Genomes:
            fout.write(f"genome {genome}\n")
            fout.write(f"trackDb {genome}/trackDb.txt\n")
        fout.write("\n")

    # make ChIP_Atlas.html
    outfile=f'{outfolder}/ChIP_Atlas.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("ChIP data from https://chip-atlas.org/ mm10 and hg38 TF\n")

    # make url.txt
    outfile=f'{outfolder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("http://upnaesrv1.epfl.ch/ChIP_Atlas/hub.txt")

    # make trackDb.txt
    for genome in Genomes:

        # get all TFs from ChIP-Atlas
        infile = f"~/ChIP_atlas/resources/experimentList_{genome}_TFs_only_QC_filtered.tab"
        Chip = pd.read_csv(infile,sep='\t',index_col=0,usecols=(0,3))
        Chip.columns=['tf']
        # sort by tf name
        Chip = Chip.sort_values('tf',axis=0)

        # make trackDb.txt
        outfile=f'{outfolder}/{genome}/trackDb.txt'
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

                # make subtracks
                i = 0
                n = sum(Chip.tf==tf)
                G = min(n+159,240)
                for idx in Chip.loc[Chip.tf==tf].index:
                    if os.path.exists(f'resources/tracks/{genome}/{idx}.bw'):
                        try:
                            pyBigWig.open(f'resources/tracks/{genome}/{idx}.bw')
                        except:
                            print('file corrupted')
                        finally:
                            i+=1

                            fout.write(f"\ttrack ChIP_{tf}_{i}\n")
                            fout.write(f"\tbigDataUrl http://upnaesrv1.epfl.ch/ChIP_Atlas/tracks/{genome}/{idx}.bw\n")
                            fout.write(f"\tshortLabel {tf}_{i}\n")
                            fout.write(f"\tlongLabel Chip_{tf}_{i}\n")
                            fout.write(f"\tparent ChIP_{tf}\n")
                            fout.write("\ttype bigWig\n")
                            fout.write(f"\tcolor 255,{G},255\n")
                            fout.write("\n")

                fout.write("###################\n\n")
