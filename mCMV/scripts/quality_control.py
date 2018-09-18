# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Initial quality control of samples.
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
def plot_cumulative(x, ax=None, **kwargs):
    '''Plot cumulative histogram'''
    if ax is None:
        fig, ax = plt.subplots()

    if isinstance(x, pd.DataFrame):
        x = x.values

    x = np.sort(x)
    y = 1.0 - np.linspace(0, 1, len(x))
    plot_kwargs = {
            'lw': 2,
            }
    plot_kwargs.update(kwargs)
    ax.plot(x, y, **plot_kwargs)

    return ax


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--save', action='store_true',
                        help='Store filtered cells dataframe of counts to file')
    parser.add_argument('--n-reads-min', type=int, default=15000,
                        help='Minimal number of reads for good cells')
    args = parser.parse_args()

    expname = 'mouse_mCMV_1'
    samplename = 'all'

    print('Read counts table')
    fn = get_count_filenames(expname, samplename, fmt='dataframe')
    counts = CountsTable(1e6 * pd.read_csv(fn, sep='\t', index_col=0))
    counts._normalized = 'counts_per_million'

    print('Read gene metadata')
    fn = '../../data/mouse_mCMV_1/feature_metadata.tsv'
    featuresheet = FeatureSheet(pd.read_csv(fn, sep='\t', index_col=0))

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

    def plot_qcs(ds, virus_threshold=60):
        fig, axs = plt.subplots(2, 3, figsize=(12, 7))
        axs = axs.ravel()

        # Number of reads
        ax = axs[0]
        col = 'n_reads'
        col_label = 'Number of reads'
        plot_cumulative(0.1 + ds.samplesheet[col], label='all cells', color='k', ax=ax)
        for sn, datum in ds.samplesheet[[col, 'biosample']].groupby('biosample'):
            x = 0.1 + datum[col]
            plot_cumulative(x, label=sn, ax=ax)
        ax.grid(True)
        ax.set_xlabel(col_label)
        ax.set_xlim(xmin=0.9 * ds.samplesheet[col].min())
        ax.set_xscale('log')

        # Number of genes
        ax = axs[1]
        col = 'n_genes_1+'
        col_label = 'Number of Genes (1+)'
        plot_cumulative(0.1 + ds.samplesheet[col], label='all cells: {:}'.format(ds.samplesheet.shape[0]), color='k', ax=ax)
        for sn, datum in ds.samplesheet[[col, 'biosample']].groupby('biosample'):
            x = 0.1 + datum[col]
            plot_cumulative(x, label='{:}: {:}'.format(sn, len(x)), ax=ax)
        ax.grid(True)
        ax.set_xlabel(col_label)
        ax.set_xlim(xmin=0.9 * ds.samplesheet[col].min())
        ax.set_xscale('log')
        ax.legend(loc='lower left', fontsize=8)

        # Number of genes
        ax = axs[2]
        col = 'n_genes_3+'
        col_label = 'Number of Genes (3+)'
        plot_cumulative(0.1 + ds.samplesheet[col], label='all cells', color='k', ax=ax)
        for sn, datum in ds.samplesheet[[col, 'biosample']].groupby('biosample'):
            x = 0.1 + datum[col]
            plot_cumulative(x, label=sn, ax=ax)
        ax.grid(True)
        ax.set_xlabel(col_label)
        ax.set_xlim(xmin=0.9 * ds.samplesheet[col].min())
        ax.set_xscale('log')

        # Number of virus reads
        ax = axs[3]
        col = 'n_reads_virus'
        col_label = 'Number of virus reads'
        plot_cumulative(0.1 + ds.samplesheet[col], label='all cells', color='k', ax=ax)
        for sn, datum in ds.samplesheet[[col, 'biosample']].groupby('biosample'):
            x = 0.1 + datum[col]
            plot_cumulative(x, label=sn, ax=ax)
        ax.plot([virus_threshold] * 2, [0, 1], lw=1.5, color='k', alpha=0.7, ls='--')
        ax.grid(True)
        ax.set_xlabel(col_label)
        ax.set_xlim(xmin=0.09)
        ax.set_xscale('log')

        # Housekeeping
        gnames = ['Actb', 'Gapdh']
        ind = ds.featuresheet.loc[ds.featuresheet['GeneName'].isin(gnames)].index
        dsind = Dataset(
                counts_table=ds.counts.loc[ind],
                samplesheet=ds.samplesheet,
                featuresheet=ds.featuresheet.loc[ind],
                )
        dsind.rename(axis='features', column='GeneName', inplace=True)
        for gname, ax in zip(gnames, axs[4:]):
            dd = dsind.counts.log().loc[[gname]].T
            dd['biosample'] = dsind.samplesheet['biosample']
            sns.violinplot(
                    data=dd,
                    x='biosample',
                    y=gname,
                    ax=ax,
                    zorder=10,
                    )
            ax.grid(True)
            ax.set_ylim(-1, 5)
            ax.set_yticks(np.arange(-1, 6))
            ax.set_yticklabels([
                    '$0$',
                    '$1$',
                    '$10$',
                    '$10^2$',
                    '$10^3$',
                    '$10^4$',
                    '$10^5$',
                    ])
            ax.set_ylabel('{:} per million reads'.format(gname))
            ax.set_xlabel('')
            for tk in ax.get_xticklabels():
                tk.set_rotation(300)

        return fig

    # QC all cells
    fig = plot_qcs(ds)
    fig.suptitle('All cells')
    plt.tight_layout(rect=(0, 0, 1, 0.96))

    # Filter cells and reQC
    n_reads_min = args.n_reads_min
    dsgood = ds.query_samples_by_metadata('n_reads > @n_reads_min', local_dict=locals())
    fig = plot_qcs(dsgood)
    fig.suptitle('Only cells with {:d}+ reads'.format(n_reads_min))
    plt.tight_layout(rect=(0, 0, 1, 0.96))

    ## Probably not worth saving, we only trash ~2,500 out of ~8,700 cells
    #if args.save:
    #    dsgood.counts.to_csv(
    #        '../../data/mouse_mCMV_1/all/dataframe_n_reads_min_{:}.tsv'.format(n_reads_min),
    #        index=True, sep='\t')

    plt.ion()
    plt.show()
