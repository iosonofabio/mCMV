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
    dsg = ds.query_features_by_metadata('(module in ("early", "mid", "late")) | (expression_geometric_mean > 100)')
    dsg.query_samples_by_metadata('moi in ("low", "high")', inplace=True)
    dsg.query_samples_by_metadata('virus_reads_per_million >= 100', inplace=True)

    print('Make metadata for module expression')
    dsg.samplesheet['exp_module_early'] = dsg.counts.loc[ds.featuresheet['module'] == 'early'].sum(axis=0)
    dsg.samplesheet['exp_module_mid'] = dsg.counts.loc[ds.featuresheet['module'] == 'mid'].sum(axis=0)
    dsg.samplesheet['exp_module_late'] = dsg.counts.loc[ds.featuresheet['module'] == 'late'].sum(axis=0)

    corrs = dsg.correlation.correlate_features_phenotypes(
        phenotypes=['exp_module_early', 'exp_module_mid', 'exp_module_late', 'virus_reads_per_million'])
    corrs = corrs.loc[dsg.featuresheet['Organism'] == 'mouse']

    print('Plot host/virus module expression for top correlates/anticorrelates')
    ctop = corrs['virus_reads_per_million'].nlargest(6).index
    cbot = corrs['virus_reads_per_million'].nsmallest(6).index

    fig, axs = plt.subplots(2, 5, figsize=(11.5, 5), sharex=True)
    x = dsg.samplesheet['log_virus_reads_per_million'].values
    for cnames, axsr, tcorner, fcolor in zip([ctop, cbot], axs, ['lower right', 'lower left'], ['tomato', 'slateblue']):
        dsg.counts.log(inplace=True)
        fits = dsg.fit.fit_single(cnames.tolist(), ['log_virus_reads_per_million'], model='threshold-linear')[:, 0]
        dsg.counts.unlog(inplace=True)
        for cname, ax, fit in zip(cnames, axsr, fits):
            y = np.log10(0.1 + dsg.counts.loc[cname]).values

            ## Fit
            #b, i, s, _ = fit.values
            #xfit = np.linspace(-2, 6, 1000)
            #yfit = i + s * xfit
            #yfit[xfit <= ((b - i) / s)] = b
            #ax.plot(xfit, yfit, ls='--', color='darkred', lw=2, zorder=11)

            # Running window
            w = 0.2
            yw = [(y[(x > (2 + i * w)) & (x <= (2 + (i + 1) * w * 2.5))]).mean() for i in range(15)]
            xw = 2 + w * (0.5 + np.arange(15))
            ax.plot(xw, yw, ls='--', color='darkred', lw=2, zorder=11)

            # KDE
            sns.kdeplot(x, y, zorder=9, cmap='magma', ax=ax, alpha=0.3)
            # scatter
            ax.scatter(x, y, s=10, zorder=10, alpha=0.05)
            for tk in [0, 1, 2, 3, 4, 5]:
                ax.plot([-2, 6], [tk] * 2, lw=0.5, color='grey', alpha=0.6, zorder=1)
                ax.plot([tk] * 2, [-2, 6], lw=0.5, color='grey', alpha=0.6, zorder=1)
            ax.set_title(cname)
            ax.add_patch(mpl.patches.Rectangle(
                (0, 0), 1, 1,
                facecolor='none',
                edgecolor=fcolor,
                lw=4,
                transform=ax.transAxes,
                ))
            ax.set_xticks(np.arange(-1, 6))
            ax.set_yticks(np.arange(-1, 6))
            ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_xlim(1.9, 5.5)
            ax.set_ylim(-1.1, y.max() * 1.1)
            label = '$\\rho = {:.2f}$'.format(
                    corrs.loc[cname, 'virus_reads_per_million'],
                    )
            if tcorner == 'lower right':
                ax.text(0.98, 0.02, label, ha='right', va='bottom',
                        transform=ax.transAxes)
            elif tcorner == 'lower left':
                ax.text(0.02, 0.02, label, ha='left', va='bottom',
                        transform=ax.transAxes)
    fig.text(0.5, 0.01, 'virus module per million RNA', ha='center', va='bottom')
    fig.text(0.01, 0.5, 'gene counts per million RNA', ha='left', va='center', rotation=90)
    plt.tight_layout(rect=(0.03, 0.03, 1, 1))

    #print('Plot distribution of virus reads')
    #fig, ax = plt.subplots()
    #x = np.sort(dsg.samplesheet['virus_reads_per_million'].values)
    #y = 1.0 - np.linspace(0, 1, len(x))
    #ax.plot(x, y, lw=2)
    #ax.set_xlabel('Virus reads per million')
    #ax.set_ylabel('Cumulative')
    #ax.grid(True)
    #plt.tight_layout()

    plt.ion()
    plt.show()

