# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/06/18
content:    Include 3'UTR annotations into GTF file.
'''
import os
import csv
import numpy as np
import pandas as pd



if __name__ == '__main__':

    print('Load 3\'UTR annotations from TSV')
    fn = '../../data/mouse_mCMV_1/mCMV_genome/mCMV_gene_annotations.tsv'
    annotsv = pd.read_csv(fn, sep='\t', index_col=0)

    print('Load GTF of the CDSes')
    fn = '../../data/mouse_mCMV_1/mCMV_genome/genes_mCMV.gtf'
    anno_cds = pd.read_csv(fn, sep='\t', header=None)

    print('Mix 3\'UTRs into GTF')
    gtf = []
    for _, line in anno_cds.iterrows():
        if line[2] not in ('gene', 'transcript'):
            gtf.append(line)
            continue
        # Positions of the genomic location
        start, end, strand = line[[3, 4, 6]]
        locstr = '{:}..{:}'.format(start, end)
        anno = annotsv.loc[annotsv['Genomic location'] == locstr]
        if anno.shape[0] != 1:
            raise ValueError('Gene not found: {:}'.format(locstr))
        anno = anno.iloc[0]

        # Reannotate gene end to include 3'UTR
        if anno["3'UTR length"] > 0:
            if strand == '+':
                line[4] = anno['End gene']
            else:
                line[3] = anno['End gene']
        gtf.append(line)

        if line[2] == 'gene':
            continue

        # Include additional exon
        if anno["3'UTR length"] > 0:
            line = line.copy()
            line[2] = 'exon'
            if strand == '+':
                line[3] = anno['Start 3\'UTR']
                line[4] = anno['End 3\'UTR (estimate)']
            else:
                line[4] = anno['Start 3\'UTR']
                line[3] = anno['End 3\'UTR (estimate)']
            line[8] += '; exon_type "3UTR (estimated)"'
            gtf.append(line)

    gtf = pd.DataFrame(gtf)
    for col in anno_cds.columns:
        gtf.loc[:, col] = gtf.loc[:, col].astype(anno_cds[col].dtype)

    print('Write new GTF to file')
    fn_out = '../../data/mouse_mCMV_1/mCMV_genome/genes_mCMV_with3UTR.gtf'
    # NOTE: to_csv quotes all text (!)
    gtf.to_csv(fn_out, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
