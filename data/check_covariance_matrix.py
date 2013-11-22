
"""
Show how the PMBEC covariance matrix was inverted.


"""
import math

import numpy as np
from scipy.stats import pearsonr

import matplotlib.pyplot as plt

def read_matrix(fname, header=False):
    mat = []
    idx = 0
    if header == True: idx=1
    lines=open(fname,'r').readlines()
    for l in lines[idx:]:
        l = l.strip()
        row = l.split()
        row = map(float, row[idx:])
        mat.append(row)
        #print len(row)
    mat = np.array(mat)
    assert mat.shape == (20,20)
    return mat


def compare_eig(ma,mb):
    """
    Compare eigenvalues and eigenvectors of the two matrices.
    """
    wa,va = np.linalg.eig(ma)
    wb,vb = np.linalg.eig(mb)

    result = sorted(zip(wa,va.transpose()), key=lambda row:row[0])
    row_list = [row[1] for row in result]
    v = np.array(row_list)
    va = v.transpose()

    result = sorted(zip(wb,vb.transpose()), key=lambda row:row[0])
    row_list = [row[1] for row in result]
    v = np.array(row_list)
    vb = v.transpose()


    #Plot eigenvalues: Exclude one with zero eigenvalue.
    plt.plot(sorted(wa)[1:], sorted(wb)[0:-1], 'ko')
    x = [min(sorted(wa)[1:]), max(sorted(wa)[1:])]
    plt.plot(x,x, '-')
    plt.show()



def analyze_pmbec(offset=0.05):
    m = read_matrix(fname_pmbec,header=True) 
    m_inv_orig = read_matrix(fname_pmbec_inv,header=False) 

    m_inv = np.linalg.inv(m)

    mb = m + offset
    #mb = m + np.eye(20)*offset  
    mb_inv = np.linalg.inv(mb)

    xref = m_inv_orig.flatten()
    xa = m_inv.flatten()
    xb = mb_inv.flatten()
  
    #Plot correlations:
    #=================
    plt.plot(xref,xb, 'ko')
    xline = [min(xb),max(xb)]
    plt.plot(xline,xline, 'r-')
    plt.xlabel('Inverse cov matrix hardcoded in the c++ code')
    plt.ylabel('Inverse cov matrix after adding 0.05 to all elements')
    plt.show()

    print 'M = PMBEC covariance matrix.'
    print 'det(M) = ', np.linalg.det(m)
    print 'det(M + %s)'%(offset,), '=', np.linalg.det(mb)

    print 'Correlation between the inverse matrix and pmbec w/ offset added:', pearsonr(xref, xb)[0]
    return m

def analyze_blosum():
    m = read_matrix(fname_blosum, header=True)
    m_inv = np.linalg.inv(m)
    print 'M = BLOSUM62 matrix'
    print 'det(M) ', np.linalg.det(m)
    return m, m_inv

def get_condition_number(m):
    """
    http://mathworld.wolfram.com/ConditionNumber.html
    """
    u,s,v = np.linalg.svd(m);
    cond_num = max(s)/min(s)
    return cond_num

if __name__ == '__main__':
    fname_pmbec = 'cov_matrix.pmbec.mat'
    fname_blosum = 'ncbi.blosum62.sij.mat'
    fname_pmbec_inv = 'InverseCovariance.txt' #This was used in the smmpmbec code.


    m_inv = read_matrix(fname_pmbec_inv,header=False)
    m = analyze_pmbec(offset=0.05)
    #m, m_inv = analyze_blosum()


    
