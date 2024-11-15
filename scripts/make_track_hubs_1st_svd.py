import pandas as pd
import os
import pyBigWig


if __name__ == '__main__':

    # Parameters
    Genomes = ['mm10','hg38']
    
    # Track hub name and url
    track_hub_name = "ChIP_Atlas"
    track_hub_url = f"https://sv-open.epfl.ch/upnae-public/sites/{track_hub_name}"

    # Create track hub folder
    track_hub_folder = f"/data/shared/sv-open/sites/{track_hub_name}"
    if not os.path.exists(track_hub_folder):
        os.makedirs(track_hub_folder)

    # make track hub file
    outfile=f'{track_hub_folder}/hub.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"hub {track_hub_name}\n")
        fout.write(f"shortLabel {track_hub_name}\n")
        fout.write("longLabel A ChIP-seq atlas of TFs in mouse and human tissues\n")
        fout.write("genomesFile genomes.txt\n")
        fout.write("email jeremie.breda@epfl.ch\n")
        fout.write(f"descriptionUrl {track_hub_name}.html\n")
        fout.write("\n")

    # make genomes.txt
    outfile=f'{track_hub_folder}/genomes.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        for genome in Genomes:
            fout.write(f"genome {genome}\n")
            fout.write(f"trackDb {genome}/trackDb.txt\n")
            fout.write("\n")

    # make html page
    outfile=f'{track_hub_folder}/{track_hub_name}.html'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write("ChIP signal and peaks from https://chip-atlas.org/ - mm10 and hg38 transcrtiption factors\n")

    # make url.txt
    outfile=f'{track_hub_folder}/url.txt'
    with open(outfile,'w', encoding="utf-8") as fout:
        fout.write(f"{track_hub_url}/hub.txt\n")

    # make trackDb.txt
    for genome in Genomes:
        print(genome)

        # get all TFs from ChIP-Atlas
        infile = f"~/ChIP_atlas/resources/experimentList_v3_{genome}_TFs_only_QC_filtered.tab"
        Chip = pd.read_csv(infile,sep='\t',index_col=0,usecols=(0,3))
        Chip.columns = ['tf']
        # sort by tf name
        Chip = Chip.sort_values('tf',axis=0)

        # make trackDb.txt
        os.makedirs(f"{track_hub_folder}/{genome}", exist_ok=True)
        outfile=f'{track_hub_folder}/{genome}/trackDb.txt'
        with open(outfile,'w', encoding="utf-8") as fout:

            # BigWig ChIP signal for each TF in ChIP-Atlas
            for tf in Chip.tf.unique():

                # BigWig tracks
                fout.write(f"track ChIP_{tf}\n")
                fout.write("type bigWig\n")
                fout.write(f"shortLabel {tf}\n")
                fout.write(f"longLabel ChIP-seq {tf}\n")
                fout.write(f"bigDataUrl {track_hub_url}/tracks/{genome}/bw/{tf}.bw\n")
                fout.write("visibility hide\n")
                fout.write("maxHeightPixels 500:100:8\n")
                fout.write("autoScale on\n")
                fout.write("\n")

                fout.write("###################\n\n")
