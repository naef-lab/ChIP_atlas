import numpy as np
import pandas as pd
import h5py
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Plot histogram of experiments QC')
    parser.add_argument('--infile_chip'
        ,required=True
        ,type=str
        ,help="TF Chip signal")
