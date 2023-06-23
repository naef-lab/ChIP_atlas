import json
import pandas as pd
import numpy as np
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--infile_complex'
        ,required=True
        ,type=str
        ,help="TF Complex table")
    parser.add_argument('--infile_chip'
        ,required=True
        ,type=str
        ,help="List of TFs in chip experiments")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="hdf5 with SVD and pearson's correlation")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_argument()

    # read tables
    complex_table = pd.read_csv(args.infile_complex,sep='\t')
    chip_tf = np.unique( np.array(pd.read_csv(args.infile_chip,sep='\t',header=None,usecols=[3])[3].values ))

    # declate matrix and complexes dict.
    interaction_matrix = np.zeros([len(chip_tf),len(chip_tf)])
    complexes = {}
    for complex_id in complex_table['Complex ac'].values:

        # Open complex json
        filein = f'resources/complex/json/{complex_id}.json'
        with open(filein) as f:
            d = json.load(f)
        members = []
        stoichiometry = []
        for i in range( len(d['data']) ):

            # Get all interactors of the complex
            if d['data'][i]['object'] == 'interactor':
                members.append( d['data'][i]['label'] )
            # get interaction stoichiometry
            elif d['data'][i]['object'] == 'interaction':
                for j in range( len(d['data'][i]['participants']) ):
                    if 'stoichiometry' in d['data'][i]['participants'][j].keys():
                        stoichiometry.append( int(d['data'][i]['participants'][j]['stoichiometry']) )
                    else:
                        stoichiometry.append(0)
            else:
                print(f"unknown object: {d['data'][i]['object']}")

        # get only proteins in my chip db
        proteins = list( set(members).intersection(set(chip_tf)) )    

        if len( proteins ) >= 1:

            # fill in complexes dict            
            complexes[complex_id] = {}
            complexes[complex_id]['members'] = members
            complexes[complex_id]['stoichiometry'] = stoichiometry
            
            # check for homodimer
            if (len(proteins) == 1) & (len(members) < len(stoichiometry)):
                # fill in interaction matrix
                ii = np.where(proteins[0]==chip_tf)
                interaction_matrix[ii,ii] += 1

            # heterodimer
            else:
                # fill in interaction matrix
                for pi in range(len(proteins)-1):
                    for pj in range(pi+1,len(proteins)):
                        ii = np.where(proteins[pi]==chip_tf)
                        jj = np.where(proteins[pj]==chip_tf)

                        interaction_matrix[ii,jj] += 1
                        interaction_matrix[jj,ii] += 1

        
    # Write matrix to output file
    pd.DataFrame(interaction_matrix,columns=chip_tf,index=chip_tf).to_csv(args.outfile,sep="\t")