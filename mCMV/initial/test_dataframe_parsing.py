# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Make dataframe of all and try parsing.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, SampleSheet, FeatureSheet


# Functions


# Script
if __name__ == '__main__':

    from mCMV.filenames import get_count_filenames
    from mCMV.io import read_dataframe
    from mCMV.counts_table import CountsTableSparse

    pa = argparse.ArgumentParser(description='Process some integers.')
    pa.add_argument('--samples', nargs='+', type=int,
                    help='What samples to use (1-8)')
    args = pa.parse_args()

    expname = 'mouse_mCMV_1'
    samplenames = [
            '1-uninfected',
            '2-uninfected',
            '3-mock',
            '4-mock',
            '5-low',
            '6-low',
            '7-high',
            '8-high',
            ]
    if args.samples is not None:
        samplenames = [samplenames[i-1] for i in args.samples]

    for sn in samplenames:
        fn = get_count_filenames(expname, sn, fmt='dataframe')
        df = read_dataframe(fn)
        print('{:}, n. genes: {:}, n. cells: {:}'.format(
            sn,
            df.shape[0],
            df.shape[1],
            ))
