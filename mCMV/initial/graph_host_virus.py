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

    print('Restrict to even fewer genes for plotting the graph')
    dsg = ds.query_features_by_metadata('(Organism == "mCMV") | (expression_geometric_mean > 100)')
    dsv = ds.query_features_by_metadata('Organism == "mCMV"')

    def plot_graph(dsg, ax=None, order=None, edges=None):
        print('Get knn of genes')
        import igraph as ig
        # Pearson's to start (simplest)
        #ds.counts.log(inplace=True)
        if edges is None:
            print('Calculate knn')
            n_knn = 5
            knn = dsg.graph.lshknn(axis='features', n_neighbors=n_knn, slice_length=6)
            # knn is symmetric, i.e. we can skip half
            edges = list(zip(knn.row[::2], knn.col[::2]))
        G = ig.Graph(edges)

        print('Calculate graph layout')
        layout_algorithm = 'drl'
        layout = G.layout(layout_algorithm)
        vs = pd.DataFrame(
                np.array(layout),
                index=dsg.featurenames,
                columns=['dim1', 'dim2'])

        print('Plot graph layout')
        if order is None:
            modules_unique = np.unique(dsg.featuresheet['module'])
        else:
            modules_unique = np.array(order)
        colors = sns.color_palette(n_colors=len(modules_unique))

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 7))

        for modname, color in zip(modules_unique, colors):
            vsi = vs.loc[dsg.featuresheet['module'] == modname]
            vsi.plot(kind='scatter',
                     x='dim1',
                     y='dim2',
                     ax=ax,
                     grid=False,
                     color=color,
                     label=modname,
                     alpha=0.4,
                     zorder=5,
                     )
        for edge in edges:
            x0, y0 = layout[edge[0]]
            x1, y1 = layout[edge[1]]
            lw = 0.5
            alpha = 0.02
            ax.plot([x0, x1], [y0, y1], lw=lw, color='k', alpha=alpha, zorder=4)

        ax.set_axis_off()
        ax.legend(loc='best', title='Module:')

        return {
            'edges': edges,
            'graph': G,
            'vs': vs,
            'ax': ax,
            'modules_unique': modules_unique,
            }

    def plot_neighbors_heatmap(dsv, d1, order=None, ax=None, normalized=True):
        if ax is None:
            fig, ax = plt.subplots()

        modules_unique = d1['modules_unique']
        modulesi = np.empty(dsv.n_features, int)
        for imu, mu in enumerate(modules_unique):
            modulesi[dsv.featuresheet['module'].values == mu] = imu
        edges = np.array(d1['edges'], int)
        ind = modulesi[edges]
        n_modules = len(modules_unique)
        n_neighbors = np.zeros((n_modules, n_modules), int)
        for i, j in ind:
            n_neighbors[i, j] += 1
            n_neighbors[j, i] += 1
        n_neighbors = pd.DataFrame(
                n_neighbors,
                index=modules_unique,
                columns=modules_unique,
                )
        frac_neighbors = (1.0 * n_neighbors / n_neighbors.sum(axis=0)).T

        # Normalize by the number of genes (those are the potential neighbors)
        if normalized:
            n_genes_per_module = np.array([(dsv.featuresheet['module'] == m).sum() for m in modules_unique])
            frac_genes_per_module = 1.0 * n_genes_per_module / n_genes_per_module.sum()
            frac_neighbors *= 1.0 / frac_genes_per_module

        if order is not None:
            frac_neighbors = frac_neighbors.loc[order]
            frac_neighbors = frac_neighbors.loc[:, order]

        sns.heatmap(
                frac_neighbors, cmap='viridis', ax=ax,
                linewidths=0.5,
                )
        ax.set_title('Fraction of neighbors')

        return ax

    ## Plot virus graph and virus-host graph with neighbours proportions
    #fig, axs = plt.subplots(
    #        2, 2, figsize=(10, 8),
    #        gridspec_kw={'height_ratios': (1.5, 1)},
    #        )
    #axs = axs.ravel()
    #d1 = plot_graph(dsv, ax=axs[0])
    #d2 = plot_graph(dsg, ax=axs[1], order=['early', 'mid', 'late', 'other_virus', 'mouse'])
    #axs[0].set_title('knn graph of virus genes')
    #axs[1].set_title('knn graph of virus and host genes')
    #plot_neighbors_heatmap(dsv, d1, ax=axs[2])
    #plot_neighbors_heatmap(dsg, d2, ax=axs[3], order=['early', 'mid', 'late', 'other_virus', 'mouse'])
    #plt.tight_layout()

    def plot_neighbors_heatmap_nonsym(dsv, d1, source_cat, dest_cats, ax=None, normalized=True):
        modules_unique = d1['modules_unique']
        modulesi = np.empty(dsv.n_features, int)
        for imu, mu in enumerate(modules_unique):
            modulesi[dsv.featuresheet['module'].values == mu] = imu
        edges = np.array(d1['edges'], int)
        ind = modulesi[edges]
        n_modules = len(modules_unique)
        n_neighbors = np.zeros((n_modules, n_modules), int)
        for i, j in ind:
            n_neighbors[i, j] += 1
            n_neighbors[j, i] += 1
        n_neighbors = pd.DataFrame(
                n_neighbors,
                index=modules_unique,
                columns=modules_unique,
                )
        n_neighbors_mat = n_neighbors.loc[dest_cats, [source_cat]]

        frac_neighbors = (1.0 * n_neighbors_mat / n_neighbors_mat.sum())

        # Normalize by the number of genes (those are the potential neighbors)
        if normalized:
            n_genes_per_module = np.array([(dsv.featuresheet['module'] == m).sum() for m in dest_cats])
            frac_genes_per_module = 1.0 * n_genes_per_module / n_genes_per_module.sum()
            frac_neighbors = (frac_neighbors.T / frac_genes_per_module).T

        sns.heatmap(
                frac_neighbors, cmap='viridis', ax=ax,
                linewidths=0.5,
                )
        ax.set_title('Fraction of neighbors')

        return ax

    def get_knn_symmetric(dsg, k=5, threshold=0.2):
        corr = dsg.correlation.correlate_features_features()
        corr_abs = np.abs(corr)
        corr_abs -= np.eye(corr_abs.shape[0])
        corr_abs[corr_abs < threshold] = 0
        edges = []
        for i in range(len(corr_abs)):
            nei = np.argpartition(corr_abs.iloc[i].values, -5)[-5:]
            nei = nei[corr_abs.iloc[i, nei] >= threshold]
            edges.extend([(i, ne) for ne in nei])
        return edges

    knn_edges_abs = get_knn_symmetric(dsg)

    print('Plot only virus-host with neighbours w/o host-host')
    dsg = ds.query_features_by_metadata('(module in ("early", "mid", "late")) | ((Organism == "mCMV") & (expression_geometric_mean > 0.05)) | (expression_geometric_mean > 100)')

    fig, axs = plt.subplots(
            1, 3, figsize=(10, 4),
            gridspec_kw={'width_ratios': [10, 10, 2]},
            )
    axs = axs.ravel()
    d2 = plot_graph(dsg, ax=axs[0], order=['early', 'mid', 'late', 'other_virus', 'mouse'], edges=knn_edges_abs)
    axs[0].set_title('knn graph of virus and host genes')
    plot_neighbors_heatmap(dsg, d2, order=['early', 'mid', 'late', 'other_virus', 'mouse'], ax=axs[1])
    plot_neighbors_heatmap_nonsym(dsg, d2, 'mouse', ['early', 'mid', 'late', 'other_virus'], ax=axs[2])

    plt.tight_layout()

    plt.ion()
    plt.show()
