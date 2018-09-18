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

from mCMV.filenames import get_count_filenames

# Functions


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--n-reads-min', type=int, default=15000,
                        help='Minimal number of reads for good cells')

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
    ind = (ds.counts > 10).sum(axis=1) >= 10
    ds.counts = ds.counts.loc[ind]

    print('Ignore genes with multiple IDs')
    from collections import Counter
    genec = Counter(ds.featuresheet['GeneName'].values)
    genes_multiple = [k for k, v in genec.items() if v > 1]
    ds.featuresheet = ds.featuresheet.loc[~ds.featuresheet['GeneName'].isin(genes_multiple)]

    print('Translate to gene names')
    ds.rename(axis='features', column='GeneName', inplace=True)

    print('Get average expression')
    ds.counts.log(inplace=True)
    ds.featuresheet.loc[:, 'expression_geometric_mean'] = ds.counts.mean(axis=1)

    print('Divide features in host and pathogen')
    features_host = ds.featurenames[ds.featuresheet['Organism'] == 'mouse']
    features_virus = ds.featurenames[ds.featuresheet['Organism'] == 'mCMV']

    print('Plot average expression of viral genes')
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(
            ds.featuresheet.loc[features_virus, 'expression_geometric_mean'].values,
            bins=np.linspace(-1, 2, 20),
            align='mid',
            )
    ax.grid(True)
    ax.set_xlabel('Average expression')
    ax.set_ylabel('# viral genes')
    ax.set_ylim(ymin=0.1)
    ax.set_yscale('log')
    ax.set_xlim(-1.1, 2.6)
    ax.set_xticks(np.arange(-1, 3))
    ax.set_xticklabels([
            '$0$',
            '$1$',
            '$10$',
            '$10^2$',
            ])
    plt.tight_layout()

    print('Check correlations between virus genes')
    dsv = ds.query_features_by_metadata('Organism == "mCMV"')
    dsv.query_samples_by_metadata('moi in ("low", "high")', inplace=True)
    dsv.query_samples_by_metadata('n_reads_virus > 0', inplace=True)
    cluster_samples = dsv.cluster.hierarchical(axis='samples', log_features=False, optimal_ordering=True)
    cluster_features = dsv.cluster.hierarchical(axis='features', log_features=False, optimal_ordering=True)
    g = dsv.plot.clustermap(
            cluster_samples=cluster_samples['linkage'],
            cluster_features=cluster_features['linkage'],
            labels_samples=False,
            annotate_samples={
                'log_virus_reads_per_million': 'viridis',
                'moi': 'Set1',
                'biosample': 'Set2',
                },
            annotate_features={
                'expression_geometric_mean': 'viridis',
                },
            colorbars=True,
            cbar_kws={'label': 'log10 expression'},
            figsize=[16.12, 9.82],
            )
    for tk in g.ax_row_colors.get_xticklabels():
        tk.set_rotation(0)
    for tk in g.ax_heatmap.get_yticklabels():
        tk.set_fontsize(4)
    g.ax_col_colors.set_yticklabels(['# virus reads', 'moi', 'biosample'])
    plt.subplots_adjust(bottom=0.03, right=0.93)
    for hli in [40.5, 56.5, 70.5, 95]:
        g.ax_heatmap.axhline(hli, lw=2, color='steelblue', zorder=20)
    for vli in [760, 1150, 1320, 1540]:
        g.ax_heatmap.axvline(vli, lw=2, color='steelblue', zorder=20)

    import json
    #modules = {
    #        '1': cluster_features['leaves'][71: 95],
    #        '2': cluster_features['leaves'][56: 71],
    #        '3': cluster_features['leaves'][40: 56],
    #        }
    #with open('../../data/mouse_mCMV_1/virus_gene_modules_1.json', 'wt') as f:
    #    json.dump(modules, f, indent=2)
    with open('../../data/mouse_mCMV_1/virus_gene_modules_1.json', 'rt') as f:
        modules = json.load(f)

    print('Plot gene expression of modules')
    virus_bins = np.array([0] + list(np.logspace(1.5, 5.5, 7)))
    virus_center = np.sqrt(np.maximum(10**0.5, virus_bins[:-1]) * virus_bins[1:])
    vrpm = dsv.samplesheet['virus_reads_per_million']
    exp = []
    frac_cells_exp = []
    for ib in range(len(virus_bins) - 1):
        if ib != len(virus_bins) - 1:
            ind = (vrpm >= virus_bins[ib]) & (vrpm < virus_bins[ib + 1])
        else:
            ind = (vrpm >= virus_bins[ib])
        sn = dsv.samplenames[ind]
        exp.append({})
        frac_cells_exp.append({})
        for modname, modgenes in modules.items():
            ex = 10**(dsv.counts.loc[modgenes, sn].values.mean() + 0.1)
            exp[-1][modname] = ex
            fc = (dsv.counts.loc[modgenes, sn].values > 0.1).mean()
            frac_cells_exp[-1][modname] = fc

    exp = pd.DataFrame(
            exp,
            index=pd.Index(virus_center, name='virus_reads_per_million'),
            )
    exp.columns.name = 'module'
    frac_cells_exp = pd.DataFrame(
            frac_cells_exp,
            index=pd.Index(virus_center, name='virus_reads_per_million'),
            )
    frac_cells_exp.columns.name = 'module'

    fig, axs = plt.subplots(2, 1, figsize=(4, 6), sharex=True)
    colors = sns.color_palette('Set1', n_colors=exp.shape[1])

    ax = axs[0]
    x = exp.index.values
    for modname, color in zip(exp.columns, colors):
        y = exp.loc[:, modname].values
        ax.plot(x, y, lw=2, color=color, label=modname)
    #ax.set_xlabel('Virus reads per million (time?)')
    ax.set_ylabel('Mean expression')
    ax.grid(True)
    ax.legend(loc='upper left', title='Module:')
    ax.set_xlim(xmin=1)
    ax.set_xscale('log')
    ax.set_ylim(ymin=0.1)
    ax.set_yscale('log')

    ax = axs[1]
    x = exp.index.values
    for modname, color in zip(exp.columns, colors):
        y = frac_cells_exp.loc[:, modname].values
        ax.plot(x, y, lw=2, color=color, label=modname)
    ax.set_xlabel('Virus reads per million (time?)')
    ax.set_ylabel('Fraction of cells expressing')
    ax.grid(True)
    ax.set_xlim(xmin=1)
    ax.set_xscale('log')

    plt.tight_layout()

    plt.ion()
    plt.show()

    sys.exit()

    print('Correlate number of virus reads with gene expression')
    corr = ds.correlation.correlate_features_features(
            features=features_host,
            features2=features_virus,
            method='spearman')

    corr_virvir = ds.correlation.correlate_features_features(
            features=features_virus,
            features2=features_virus,
            method='spearman')


    print('Plot top correlates')
    # Partial sort of the correlations
    tmpi = []
    tmpj = []
    # Positive correlations
    tmp = (corr).values.ravel()
    tmparg = tmp.argpartition(-16)[-16:][::-1]
    tmpj.extend(list(tmparg % corr.shape[1]))
    tmpi.extend(list(tmparg // corr.shape[1]))
    # Negative correlations
    tmp = (-corr).values.ravel()
    tmparg = tmp.argpartition(-16)[-16:][::-1]
    tmpj.extend(list(tmparg % corr.shape[1]))
    tmpi.extend(list(tmparg // corr.shape[1]))

    baseline = 1e6 / ds.samplesheet['n_reads']
    baseline_avg = np.log10(0.1 + baseline).mean()
    fig, axs = plt.subplots(4, 8, figsize=(23, 13), sharex=True, sharey=True)
    axs = axs.ravel()
    for ax, ih, iv in zip(axs, tmpi, tmpj):
        geneh = corr.index[ih]
        genev = corr.columns[iv]
        rho = corr.loc[geneh, genev]
        x = ds.counts.loc[genev].values
        y = ds.counts.loc[geneh].values
        avg = ds.featuresheet.loc[geneh, 'expression_geometric_mean']
        ax.scatter(x, y, s=15, alpha=0.05, label='$exp = {:.1f}$\n$\\rho = {:.2f}$'.format(avg, rho))
        sns.kdeplot(data=x, data2=y, legend=False, cmap='viridis', ax=ax)
        ax.axhline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.axvline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.grid(True)
        ax.set_xlabel(genev)
        ax.set_ylabel(geneh)
        ax.set_xlim(-1.1, 4.1)
        ax.set_xticks(np.arange(-1, 5))
        ax.set_xticklabels([
                '$0$',
                '$1$',
                '$10$',
                '$10^2$',
                '$10^3$',
                '$10^4$',
                ])
        ax.set_ylim(-1.1, 5.1)
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
        ax.legend(loc='best', framealpha=1)

    fig.text(0.01, 0.25, '$\\rho \ll 0$', rotation=90, va='center')
    fig.text(0.01, 0.75, '$\\rho \gg 0$', rotation=90, va='center')
    plt.tight_layout(rect=(0.02, 0, 1, 1))

    plt.ion()
    plt.show()
