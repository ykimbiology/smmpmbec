

"""
From raw binding affinity data for 24 MHC molecules, generate covariance matrix for 20 residue types. 
"""

import math

import numpy as np


fname_binding = '../../Background/FormattedCombinatorials_log_ic50.txt'

def read_table():
    d = []
    lines=open(fname_binding,'r').readlines()
    for l in lines[1:]:
        row = l.strip().split('\t')
        row = row[2:]  # Remove first 2 columns.
        #print row
        #print len(row)
        assert len(row) == 20
        row = map(float, row)  # For each position, residue contributions.
        row = np.array(row)
        row = row - np.mean(row)  # Centering by its mean.
        d.append(row)
    d = np.array(d)
    return d


d = read_table()

#PMBEC covariance matrix: 20x20
dcov = np.cov(d.transpose())

for row in dcov:
    print '\t'.join(map(str, row))
