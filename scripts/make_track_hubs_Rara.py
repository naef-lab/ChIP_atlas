import pandas as pd
import os
import pyBigWig


if __name__ == '__main__':

    Genomes = ['mm10']

    track_hub_folder = '/data/web/sites/ChIP_Atlas_Rara'
    track_hub_url = 'http://upnaesrv1.epfl.ch/ChIP_Atlas_Rara'
    
    #link to bw tracks

    # link to TF peaks
    for genome in Genomes:
        link_to_data = f'{track_hub_folder}/{genome}/bw'
        data_folder = f'/bigdata/jbreda/ChIP/resources/tracks/{genome}'

        if not os.path.exists(link_to_data):
            os.system(f'ln -s {data_folder} {link_to_data}')

        link_to_data = f'{track_hub_folder}/{genome}/bb'
        data_folder = f'/bigdata/jbreda/ChIP/results/{genome}/TF_peaks'
        if not os.path.exists(link_to_data):
            os.system(f'ln -s {data_folder} {link_to_data}')

    # make hub.txt
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("hub ChIP_Atlas\n")
        fout.write("shortLabel ChIP_atlas\n")
        fout.write("longLabel A ChIP-seq atlas of TFs in mouse and human tissues\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write("descriptionUrl ChIP_Altas.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{track_hub_folder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        for genome in Genomes:
            fout.write(f"genome {genome}\n")
            fout.write(f"trackDb {genome}/trackDb.txt\n")
            fout.write("\n")

    # make html page
    outfile=f'{track_hub_folder}/ChIP_Atlas.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("ChIP signal and peaks from https://chip-atlas.org/ - mm10 and hg38 transcrtiption factors\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("http://upnaesrv1.epfl.ch/ChIP_Atlas/hub.txt\n")

    # make trackDb.txt
    for genome in Genomes:
        print(genome)

        # get all TFs from ChIP-Atlas
        infile = f"~/ChIP_atlas/resources/experimentList_{genome}_TFs_only_QC_filtered.tab"
        Chip = pd.read_csv(infile,sep='\t',index_col=0,usecols=(0,3))
        Chip.columns = ['tf']
        # sort by tf name
        Chip = Chip.sort_values('tf',axis=0)

        # make trackDb.txt
        os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
        outfile=f'{track_hub_folder}/{genome}/trackDb.txt'
        with open(outfile,'w', encoding="utf-8") as fout:

            
            # BigBed ChIP Peaks for each TF in ChIP-Atlas
            for tf in Chip.tf.unique():
                fout.write(f"track Peaks_{tf}\n")
                fout.write("type bigBed 9\n")
                fout.write("itemRgb off\n")
                fout.write("spectrum on\n")
                fout.write(f"shortLabel {tf}_peak\n")
                fout.write(f"longLabel ChIP-seq peaks {tf}\n")
                fout.write(f"bigDataUrl {track_hub_url}/{genome}/bb/{tf}.bb\n")
                fout.write("visibility hide\n")
                fout.write("\n")

            # BigWig ChIP signal for each TF in ChIP-Atlas
            for tf in Chip.tf.unique():

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

                    i+=1

                    fout.write(f"\ttrack ChIP_{tf}_{i}\n")
                    fout.write(f"\tbigDataUrl {track_hub_url}/{genome}/bw/{idx}.bw\n")
                    fout.write(f"\tshortLabel {tf}_{i}\n")
                    fout.write(f"\tlongLabel Chip_{tf}_{i}\n")
                    fout.write(f"\tparent ChIP_{tf}\n")
                    fout.write("\ttype bigWig\n")
                    fout.write(f"\tcolor 255,{G},255\n")
                    fout.write("\n")

                fout.write("###################\n\n")
