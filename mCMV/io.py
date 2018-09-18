# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    IO functions.
'''
# Modules
import numpy as np
import pandas as pd



# Functions
def read_10X(fdict, make_dense=False):
    import scipy.io
    mat = scipy.io.mmread(fdict['matrix'])
    genes = pd.read_csv(
            fdict['genes'],
            sep='\t',
            squeeze=True,
            index_col=0,
            header=None,
            ).index
    barcodes = pd.read_csv(
            fdict['barcodes'],
            sep='\t',
            squeeze=True,
            index_col=0,
            header=None,
            ).index

    if not make_dense:
        return {
                'matrix': mat,
                'featurenames': genes,
                'samplenames': barcodes,
                }
    else:
        df = pd.DataFrame(
                data=mat.todense(),
                index=genes,
                columns=barcodes,
                )
        return df


def read_dataframe(fn):
    return pd.read_csv(fn, sep='\t', index_col=0)
