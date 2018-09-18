# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Filename factories.
'''
# Modules
import os


# Globals
VISCAMY_ROOT_FDN = os.getenv(
        'VISCAMY_ROOT_FDN',
        '/home/fabio/university/postdoc/mCMV/analysis/',
        )


# Functions
def get_count_filenames(expname, samplename, fmt='dataframe'):
    fdn = VISCAMY_ROOT_FDN+'data/{:}/{:}/'.format(
            expname,
            samplename,
            )

    if fmt == 'dataframe':
        return fdn+'dataframe.tsv'
    elif fmt == '10X':
        return {
                'matrix': fdn+'matrix.mtx',
                'barcodes': fdn+'barcodes.tsv',
                'genes': fdn+'genes.tsv',
                }
    else:
        raise ValueError('Format not supported: {:}'.format(fmt))
