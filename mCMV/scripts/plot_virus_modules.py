# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/06/18
content:    Study virus gene modules
'''
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

from viscamy.filenames import get_count_filenames




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

    print('Restrict to virus genes')
    dsv = ds.query_features_by_metadata('Organism == "mCMV"')
    dsv.query_samples_by_metadata("moi in ('low', 'high')", inplace=True)

    print('Log counts')
    dsv.counts.log(inplace=True)

    print('Plot double hierarchical clustering')
    g = dsv.plot.clustermap(
        cluster_samples=True,
        cluster_features=True,
        labels_samples=False,
        annotate_samples={
            'moi': 'Set1',
            },
        figsize=(12, 13),
        #colorbars=True,
        )
    g.ax_heatmap.set_xticklabels([])
    g.ax_row_dendrogram.remove()
    for tk in g.ax_heatmap.get_yticklabels():
        tk.set_fontsize(4)
    plt.subplots_adjust(-0.15, 0.03, 0.9, 0.98)

    plt.ion()
    plt.show()

