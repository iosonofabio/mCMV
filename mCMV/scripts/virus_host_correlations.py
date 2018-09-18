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

    #print('Read counts table')
    #fn = get_count_filenames(expname, samplename, fmt='dataframe')
    #counts = CountsTable(1e6 * pd.read_csv(
    #    fn,
    #    sep='\t',
    #    index_col=0,
    #    dtype={0: str},
    #    ))
    #counts._normalized = 'counts_per_million'

    #print('Read gene metadata')
    #fn = '../../data/mouse_mCMV_1/feature_metadata.tsv'
    #featuresheet = FeatureSheet(pd.read_csv(
    #    fn,
    #    sep='\t',
    #    index_col=0,
    #    dtype={0: str},
    #    ))

    #print('Read sample metadata')
    #fn = '../../data/mouse_mCMV_1/all/samplesheet.tsv'
    #samplesheet = SampleSheet(pd.read_csv(
    #    fn,
    #    sep='\t',
    #    index_col=0,
    #    dtype={0: str},
    #    ))

    #print('Build dataset')
    #ds = Dataset(
    #        counts_table=counts,
    #        featuresheet=featuresheet,
    #        samplesheet=samplesheet,
    #        )

    print('Load dataset')
    ds = Dataset(
            counts_table='combined',
            featuresheet='combined',
            samplesheet='combined',
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

    # END OF BOILERPLATE

    # Average expression
    ds.featuresheet['meanExp'] = ds.counts.mean(axis=1)

    print('Correlate number of virus reads with gene expression')
    corr = ds.correlation.correlate_features_phenotypes(phenotypes=['virus_reads_per_million'], method='spearman')
    corr.rename(columns={'virus_reads_per_million': 'rho'}, inplace=True)
    ds.featuresheet.loc[:, 'rho'] = corr.loc[:, 'rho']
    for ng in (20, 50, 100):
        top = ds.featuresheet.query('Organism == "mouse"').nlargest(ng, columns=['rho']).index
        with open('../../data/mouse_mCMV_1/results/host_correlations_top{:}.txt'.format(ng), 'wt') as f:
            f.write('\n'.join(top))
        bot = ds.featuresheet.query('Organism == "mouse"').nsmallest(ng, columns=['rho']).index
        with open('../../data/mouse_mCMV_1/results/host_correlations_bottom{:}.txt'.format(ng), 'wt') as f:
            f.write('\n'.join(bot))
    print(' '.join(top[:10]))
    print(' '.join(bot[:10]))

    print('Correlate only in infected cultures')
    dsi = ds.query_samples_by_metadata('moi in ("low", "high")', inplace=False)
    dsi.featuresheet['meanExp'] = dsi.counts.mean(axis=1)
    corri = dsi.correlation.correlate_features_phenotypes(phenotypes=['virus_reads_per_million'], method='spearman')
    corri.rename(columns={'virus_reads_per_million': 'rho'}, inplace=True)
    dsi.featuresheet.loc[:, 'rho'] = corri.loc[:, 'rho']
    for ng in (20, 50, 100):
        top = dsi.featuresheet.query('Organism == "mouse"').nlargest(ng, columns=['rho']).index
        with open('../../data/mouse_mCMV_1/results/host_correlations_infected_top{:}.txt'.format(ng), 'wt') as f:
            f.write('\n'.join(top))
        bot = dsi.featuresheet.query('Organism == "mouse"').nsmallest(ng, columns=['rho']).index
        with open('../../data/mouse_mCMV_1/results/host_correlations_infected_bottom{:}.txt'.format(ng), 'wt') as f:
            f.write('\n'.join(bot))
    print(' '.join(top[:10]))
    print(' '.join(bot[:10]))

    print('Check cell cycle')
    genes_cc = [
            'Ccna2', 'Ccnb1', 'Ccnb2', 'Ccnc', 'Ccnd1', 'Ccnd2', 'Ccnd3',
            'Ccne1', 'Ccne2', 'Ccnf', 'Ccng1', 'Ccng2', 'Ccnh', 'Ccni',
            'Ccnj', 'Ccnk', 'Ccnl1', 'Ccnl2', 'Ccno', 'Ccnt1', 'Ccnt2', 'Ccny',
            'Cdk1', 'Cdk2', 'Cdk4', 'Cdk6',
            ]
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    colors = ['steelblue' if x[1] == 'c' else 'darkred' for x in genes_cc]
    dsi.featuresheet.loc[genes_cc, ['meanExp', 'rho']].plot(
            x='meanExp',
            y='rho',
            kind='scatter',
            color=colors,
            ax=ax,
            s=20,
            )
    for gene in genes_cc:
        x = dsi.featuresheet.loc[gene, 'meanExp']
        y = dsi.featuresheet.loc[gene, 'rho']
        ax.text(x * 1.1, y + 0.02, gene)
    ax.set_xlim(0.3, 3000)
    ax.set_xscale('log')
    ymax = (np.abs(dsi.featuresheet.loc[genes_cc, 'rho'].values)).max() * 1.1
    ax.set_ylim(-ymax, ymax)
    ax.grid(True)
    plt.tight_layout()
    plt.ion()
    plt.show()


    print('Select mre cell cycle genes for dimensionality reduction')
    genes_cc2 = genes_cc + [
        'Rfc4', 'Mcm6', 'Mcm5', 'Ube2c', 'Cenpf', 'Tpx2', 'Aurka', 'Cenpe',
        'Plk1', 'Mapk13', 'Cdca8', 'Hjurp', 'Kpna2', 'Kif23', 'Cks2', 'Dtl',
        'Top2a', 'Bub1', 'Arl6ip1', 'Dlgap5', 'Ube2s', 'Nuf2', 'Hmmr', 'Cdc20',
        'Birc5',
        ]
    dscc = ds.query_features_by_name(genes_cc2, inplace=False)

    sys.exit()

    print('Plot a few notable host genes')
    from sklearn.neighbors.kde import KernelDensity
    genes = ['Uba52', 'Tmsb10', 'Lgals1', 'Sec61g', 'Ccl2', 'Cxcl1', 'Cxcl5', 'Ccl5']
    fig, axs = plt.subplots(2, 4, figsize=(11, 5), sharex=True, sharey=True)
    axs = axs.ravel()
    x = 0.1 + ds.samplesheet['virus_reads_per_million'].values
    for gene, ax in zip(genes, axs):
        y = 0.1 + ds.counts.loc[gene].values
        ax.set_title(gene)
        ax.scatter(np.log10(x), np.log10(y), s=15, alpha=0.3)
        # Plot histograms
        xband = .3
        xlefts = np.logspace(-1, 6, 14)
        maxs = []
        for xl in xlefts:
            xr = xl * (1 + xband)
            yi = np.log10(y[(x >= xl) & (x < xr)])[:, np.newaxis]
            if len(yi) < 10:
                continue

            yplot = np.log10(np.logspace(-1, 6, 1000))[:, np.newaxis]

            kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(yi)
            log_dens = np.exp(kde.score_samples(yplot))
            log_dens /= log_dens.max()

            maxs.append([np.log10(xl) + 0.5 * xband, yplot[np.argmax(log_dens), 0]])

            ax.plot(np.log10(xl) + log_dens, yplot[:, 0], lw=2)
            ax.grid(True)

        xmaxs, ymaxs = np.array(maxs).T
        ax.plot(xmaxs, ymaxs, lw=2, ls='--', color='black', alpha=0.5)

    ax.set_xlim(-1.1, 6)
    ax.set_ylim(-1.1, 6)
    fig.text(0.5, 0.02, 'virus cpm')
    axs[0].set_ylabel('gene cpm')
    plt.tight_layout(rect=(0, 0.03, 1, 1))

    plt.ion()
    plt.show()

    print('Get average expression')
    ds.counts.log(inplace=True)
    ds.featuresheet.loc[:, 'expression_geometric_mean'] = ds.counts.mean(axis=1)

    print('Plot some cases')
    samples_infected = ("5-low", "6-low", "7-high", "8-high")
    dsp = ds.query_samples_by_metadata('biosample in @samples_infected', local_dict=locals())
    dsp.featuresheet.loc[:, 'rho'] = dsp.correlation.correlate_features_phenotypes(phenotypes=['virus_reads_per_million'], method='spearman').iloc[:, 0]
    fig, axs = plt.subplots(2, 10, figsize=(23, 6), sharex=True, sharey=True)
    # Get 1 / n_reads baseline
    baseline = 1e6 / dsp.samplesheet['n_reads']
    baseline_avg = np.log10(0.1 + baseline).mean()
    baseline_std = np.log10(0.1 + baseline).std()
    group = dsp.featuresheet.groupby('Organism').get_group('mouse')['rho']
    for ax, (gname, rho) in zip(axs[0], group.nlargest(10).items()):
        x = dsp.samplesheet['log_virus_reads_per_million'].values
        y = dsp.counts.loc[gname].values
        avg = dsp.featuresheet.loc[gname, 'expression_geometric_mean']
        ax.scatter(x, y, s=15, alpha=0.05, label='$exp = {:.1f}$\n$\\rho = {:.2f}$'.format(avg, rho))
        sns.kdeplot(data=x, data2=y, legend=False, cmap='viridis', ax=ax)
        ax.axhline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.axvline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.grid(True)
        ax.set_title(gname)
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
    for ax, (gname, rho) in zip(axs[1], group.nsmallest(10).items()):
        x = dsp.samplesheet['log_virus_reads_per_million'].values
        y = dsp.counts.loc[gname].values
        avg = dsp.featuresheet.loc[gname, 'expression_geometric_mean']
        ax.scatter(x, y, s=15, alpha=0.05, label='$exp = {:.1f}$\n$\\rho = {:.2f}$'.format(avg, rho))
        sns.kdeplot(data=x, data2=y, legend=False, cmap='viridis', ax=ax)
        ax.axhline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.axvline(baseline_avg, lw=2, color='darkred', zorder=5)
        ax.grid(True)
        ax.set_title(gname)
        ax.set_xlim(-1.1, 5.5)
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

    fig.text(0.5, 0.02, 'Virus reads per million')
    fig.text(0.01, 0.5, 'Gene reads per million', rotation=90, va='center')
    fig.text(0.01, 0.25, '$\\rho \ll 0$', rotation=90, va='center')
    fig.text(0.01, 0.75, '$\\rho \gg 0$', rotation=90, va='center')
    plt.tight_layout(rect=(0.02, 0.03, 1, 1))

    plt.ion()
    plt.show()
