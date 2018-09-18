# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Make dataframe of all and try parsing.
'''
# Modules
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, SampleSheet, FeatureSheet


# Functions


# Script
if __name__ == '__main__':

    from viscamy.filenames import get_count_filenames
    from viscamy.io import read_dataframe

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

    print('Get counts from single samples')
    data = []
    cells = []
    n_cells = 0
    for sn in samplenames:
        print(sn)
        fn = get_count_filenames(expname, sn, fmt='dataframe')
        df = read_dataframe(fn)
        n_cells += df.shape[1]
        n_genes = df.shape[0]
        genes = df.index
        cells.extend([sn+'_'+bc for bc in df.columns])
        data.append(df.values)
        print('{:}, n. genes: {:}, n. cells: {:}'.format(
            sn,
            df.shape[0],
            df.shape[1],
            ))

    print('Merge into a single matrix')
    counts = np.zeros((n_genes, n_cells), float)
    i = 0
    samplesheet = []
    for datum in data:
        nc = datum.shape[1]
        counts[:, i: i + nc] = datum
        i += nc

    print('Get feature sheet')
    fn = '../../data/mouse_mCMV_1/feature_metadata.tsv'
    featuresheet = pd.read_csv(fn, sep='\t', index_col=0)

    print('Generate samplesheet')
    samplesheet = pd.DataFrame(
            data=[],
            index=cells,
            columns=[],
            )
    samplesheet['biosample'] = [sn.split('_')[0] for sn in samplesheet.index]
    samplesheet['moi'] = [bs.split('-')[1] for bs in samplesheet['biosample'].values]
    samplesheet['n_reads'] = counts.sum(axis=0).astype(int)
    samplesheet['n_reads_host'] = counts[featuresheet['Organism'] == 'mouse'].sum(axis=0).astype(int)
    samplesheet['n_reads_virus'] = counts[featuresheet['Organism'] == 'mCMV'].sum(axis=0).astype(int)
    samplesheet['n_genes_1+'] = (counts >= 1).sum(axis=0)
    samplesheet['n_genes_3+'] = (counts >= 3).sum(axis=0)
    samplesheet['n_genes_10+'] = (counts >= 10).sum(axis=0)

    print('Perform L1 normalization on counts')
    print('Number of virus reads is < 3% so we can include them')
    counts /= counts.sum(axis=0)

    print('Make counts dataframe')
    df = pd.DataFrame(
            data=counts,
            index=genes,
            columns=cells,
            )

    print('Save to file')
    df.to_csv('../../data/mouse_mCMV_1/all/dataframe.tsv', index=True, sep='\t')

