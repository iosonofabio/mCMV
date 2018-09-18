# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/03/18
content:    Try and represent host and virus genes in a knn graph.
'''
# Modules
import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, SampleSheet, FeatureSheet

from mCMV.filenames import get_count_filenames




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--n-reads-min', type=int, default=15000,
                        help='Minimal number of reads for good cells')
    parser.add_argument('--n_cpm_min_genes', nargs=2, default=(10, 10),
                        help='Only keep genes with >= x cpm in y cells.')

    args = parser.parse_args()

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

    print('Add normalized virus counts')
    ds.samplesheet['virus_reads_per_million'] = 1e6 * ds.samplesheet['n_reads_virus'] / ds.samplesheet['n_reads']
    ds.samplesheet['log_virus_reads_per_million'] = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])

    print('Filter low-quality cells')
    n_reads_min = args.n_reads_min
    ds.query_samples_by_metadata('n_reads > @n_reads_min', local_dict=locals(), inplace=True)

    print('Limit to decently expressed genes')
    ind = (ds.counts > args.n_cpm_min_genes[0]).sum(axis=1) >= args.n_cpm_min_genes[1]
    ds.counts = ds.counts.loc[ind]

    print('Ignore genes with multiple IDs')
    from collections import Counter
    genec = Counter(ds.featuresheet['GeneName'].values)
    genes_multiple = [k for k, v in genec.items() if v > 1]
    ds.featuresheet = ds.featuresheet.loc[~ds.featuresheet['GeneName'].isin(genes_multiple)]

    print('Translate to gene names')
    ds.rename(axis='features', column='GeneName', inplace=True)

    print('Add module info')
    with open('../../data/mouse_mCMV_1/virus_gene_modules_1.json', 'rt') as f:
        modules = json.load(f)
    ds.featuresheet['module'] = ds.featuresheet['Organism']
    ds.featuresheet.loc[ds.featuresheet['Organism'] == 'mCMV', 'module'] = 'other_virus'
    ds.featuresheet.loc[modules['1'], 'module'] = 'early'
    ds.featuresheet.loc[modules['2'], 'module'] = 'mid'
    ds.featuresheet.loc[modules['3'], 'module'] = 'late'

    print('Get average expression')
    ds.featuresheet.loc[:, 'expression_geometric_mean'] = 10**((np.log10(0.1 + ds.counts)).mean(axis=1)) - 0.1

    print('Restrict to genes within the modules')
    dsg = ds.query_features_by_metadata('(module in ("early", "mid", "late")) | ((Organism == "mCMV") & (expression_geometric_mean > 0.05)) | (expression_geometric_mean > 100)')

    print('Export to TSV')
    fdn = '../../data/mouse_mCMV_1/to_Paolo/'
    co = (1e-6 * dsg.counts * samplesheet['n_reads']).fillna(0).astype(int)
    co.to_csv(
            fdn+'counts_not_normalized_filtered.tsv',
            sep='\t',
            index=True,
            )
    dsg.samplesheet.to_csv(
            fdn+'samplesheet_filtered.tsv',
            sep='\t',
            index=True,
            )
    dsg.featuresheet.to_csv(
            fdn+'featuresheet_filtered.tsv',
            sep='\t',
            index=True,
            )
