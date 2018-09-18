# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/03/18
content:    Data structures for sparse count tables.
'''
# Modules
import numpy as np
import pandas as pd


# Functions
class CountsTableSparse:
    def __init__(self, matrix, featurenames, samplenames):
        self.matrix = matrix
        self.featurenames = featurenames
        self.samplenames = samplenames

    @classmethod
    def from_10X_filenames(cls, fdict):
        from mCMV.io import read_10X
        fdict = read_10X(fdict, make_dense=False)
        return cls(**fdict)

    @classmethod
    def from_10X(cls, expname, samplename):
        from mCMV.filenames import get_count_filenames
        fdict = get_count_filenames(expname, samplename, fmt='10X')
        return cls.from_10X_filenames(fdict)

    def todense(self):
        df = pd.DataFrame(
                data=self.matrix.todense(),
                index=self.featurenames,
                columns=self.samplenames,
                )
        return df
