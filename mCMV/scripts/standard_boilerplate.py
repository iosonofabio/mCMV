# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Boilerplate code to start pretty much any analysis on mCMV.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, SampleSheet, FeatureSheet

from viscamy.filenames import get_count_filenames

# Functions


# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Process some integers.')
    pa.add_argument('--save', action='store_true',
                    help='Store filtered cells dataframe of counts to file')
    pa.add_argument('--n-reads-min', type=int, default=15000,
                    help='Minimal number of reads for good cells')
    pa.add_argument('--keep2', action='store_true',
                    help='Keep sample 2-uninfected despite low quality')
    args = pa.parse_args()

    expname = 'mouse_mCMV_1'
    samplename = 'all'

    print('Read counts table')
    fn = get_count_filenames(expname, samplename, fmt='dataframe')
    counts = CountsTable(1e6 * pd.read_csv(
        fn,
        sep='\t',
        index_col=0,
        dtype={0: str},
        ))
    counts._normalized = 'counts_per_million'

    print('Read gene metadata')
    fn = '../../data/mouse_mCMV_1/feature_metadata.tsv'
    featuresheet = FeatureSheet(pd.read_csv(
        fn,
        sep='\t',
        index_col=0,
        dtype={0: str},
        ))

    print('Read sample metadata')
    fn = '../../data/mouse_mCMV_1/all/samplesheet.tsv'
    samplesheet = SampleSheet(pd.read_csv(
        fn,
        sep='\t',
        index_col=0,
        dtype={0: str},
        ))

    print('Build dataset')
    ds = Dataset(
            counts_table=counts,
            featuresheet=featuresheet,
            samplesheet=samplesheet,
            )

    print('Filter low-quality cells')
    n_reads_min = args.n_reads_min
    ds.query_samples_by_metadata('n_reads > @n_reads_min', local_dict=locals(), inplace=True)

    if not args.keep2:
        print('Filter out sample 2-uninfected (low quality)')
        ds.query_samples_by_metadata('biosample != "2-uninfected"', inplace=True)

    print('Add normalized virus counts')
    ds.samplesheet['virus_reads_per_million'] = 1e6 * ds.samplesheet['n_reads_virus'] / ds.samplesheet['n_reads']
    ds.samplesheet['log_virus_reads_per_million'] = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])

    print('Limit to decently expressed genes')
    ind = (ds.counts >= 10).sum(axis=1) >= 10
    ds.featuresheet['detected_1010'] = ind.values

    print('Ignore genes with multiple IDs')
    from collections import Counter
    genec = Counter(ds.featuresheet['GeneName'].values)
    ind = [genec[gn] for gn in ds.featuresheet['GeneName'].values]
    ds.featuresheet['nGenes'] = ind

    ds.query_features_by_metadata('((Organism == "mCMV") | detected_1010) & (nGenes == 1)', inplace=True)

    print('Translate to gene names')
    ds.featuresheet['EnsemblID'] = ds.featuresheet.index
    ds.rename(axis='features', column='GeneName', inplace=True)
