# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/06/18
content:    Recount virus reads with htseq
'''
import os
import shutil
import argparse
import numpy as np
import pandas as pd
import subprocess as sp


if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Recount mCMV reads')
    pa.add_argument('--sample', required=True,
                    help='What sample to analyze')
    pa.add_argument('--add-missing', action='store_true',
                    help='Add missing cell barcodes to the count matrix')
    args = pa.parse_args()

    fn = '../../data/mouse_mCMV_1/{:}/possorted_mapped_to_mCMV_counted.sam'.format(args.sample)
    if not os.path.isfile(fn):
        print('Count with htseq-count using SAM output')
        sp.run(
            'python ~/university/postdoc/htseq/htseq/python3/HTSeq/scripts/count.py'
            ' -f bam'
            ' --secondary-alignments ignore --supplementary-alignments ignore'
            ' -m intersection-nonempty'
            ' -s yes'
            ' -o ../../data/mouse_mCMV_1/{0}/possorted_mapped_to_mCMV_counted.sam'
            ' -r pos'
            ' ../../data/mouse_mCMV_1/{0}/possorted_mapped_to_mCMV.bam'
            ' ../../data/mouse_mCMV_1/mCMV_genome/genes_mCMV_with3UTR.gtf'.format(args.sample),
            shell=True,
            check=True,
            )

    print('Load gene ids')
    fn = '../../data/mouse_mCMV_1/{:}/genes.tsv'.format(args.sample)
    featuresheet = pd.read_csv(fn, sep='\t', index_col=0)
    virus_genes = featuresheet.index[~featuresheet.index.str.startswith('ENS')]
    virus_genes = virus_genes.append(pd.Index(
        ['__ambiguous_virus',
         '__no_feature',
         '__alignment_not_unique',
         ]))

    print('Load cell names')
    fn = '../../data/mouse_mCMV_1/{:}/barcodes.tsv'.format(args.sample)
    samplenames = pd.read_csv(fn, sep='\t', squeeze=True, header=None).values

    # Get cell barcode and assignment from htseq SAM output
    counts = pd.DataFrame(
            data=np.zeros((len(virus_genes), len(samplenames))),
            index=virus_genes,
            columns=samplenames,
            )
    print('Parsing htseq-count output SAM file')
    fn = '../../data/mouse_mCMV_1/{:}/possorted_mapped_to_mCMV_counted.sam'.format(args.sample)
    with open(fn, 'rt') as samfile:
        for il, line in enumerate(samfile):
            if ((il + 1) % 50000) == 0:
                print(il+1, 'reads parsed')
            tags = line.rstrip('\n').split('\t')[11:]
            cb_found = False
            as_found = False
            for tag in tags:
                if tag.startswith('CR:Z:'):
                    cell_barcode = tag[5:]+'-1'
                    cb_found = True
                elif tag.startswith('XF:Z:'):
                    feature = tag[5:]
                    as_found = True
                if cb_found and as_found:
                    break

            # Add new barcodes to the matrix for now
            if cell_barcode not in counts.columns:
                if args.add_missing:
                    counts[cell_barcode] = 0
                else:
                    continue

            if 'ambiguous' in feature:
                counts.loc['__ambiguous_virus', cell_barcode] += 1
                continue

            counts.loc[feature, cell_barcode] += 1
        print(il+1, 'reads parsed')

    print('Load the other counts to merge')
    fn = '../../data/mouse_mCMV_1/{:}/dataframe.tsv'.format(args.sample)
    counts_all = pd.read_csv(fn, sep='\t', index_col=0)

    fnbak = '../../data/mouse_mCMV_1/{:}/dataframe.tsv.bak'.format(args.sample)
    if not os.path.isfile(fnbak):
        print('Make backup')
        shutil.copy(fn, fnbak)

    print('Merge counts')
    counts = counts.loc[:, counts_all.columns]
    counts = counts.loc[[i for i in counts_all.index if i in counts.index[:-3]]]
    co = counts_all.values
    co[-counts.shape[0]:] = counts.values
    counts_all = pd.DataFrame(
            data=co.astype(int),
            index=counts_all.index,
            columns=counts_all.columns,
            )
    counts_all.to_csv(fn, sep='\t', index=True)
