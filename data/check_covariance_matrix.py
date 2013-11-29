
"""
Show how the PMBEC covariance matrix was inverted.


"""
import math
import random

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

def write_matrix(m, fname):
    """
    """
    aalist = list('ACDEFGHIKLMNPQRSTVWY')
    f=open(fname,'w')
    f.write('* '+' '.join(aalist)+'\n')
    for aa,row in zip(aalist,m):
        row = [aa]+list(row)
        line = ' '.join(map(str, row))
        f.write(line+'\n')
    f.close()


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
    m = (m + m.transpose())/2.0  # To ensure that the matrix is symmetric. This step does not change results.
    m_inv = np.linalg.inv(m)

    w,v = np.linalg.eig(m)
    #w.sort()

    #Numerically determine whether x*M*x >= for randomly sampled vectors, x:
    n = 10000
    xlist = []
    for i in range(n):
        x = [random.random() for j in range(20)]
        x = np.array(x)
        x = x - np.mean(x) # x = centered vector of length 20.
        x.shape = (1,20)
        xvalue = x.dot(m).dot(x.transpose())
#        xvalue = x.dot(m_inv).dot(x.transpose())
        xlist.append(xvalue[0][0])
    plt.hist(xlist); plt.xlabel('x*M*x'); plt.ylabel('freq'); plt.show()
    print 'Properties of the blosum matrix'
    print '======================================'
    print 'M = BLOSUM62 matrix'
    print 'det(M) = ', np.linalg.det(m)
    print 'Sorted_eigenvalues = ', w
    print '======================================'
    print '1) The matrix is invertible because its determinant is not zero. See: [http://en.wikipedia.org/wiki/Invertible_matrix] '
    print 'This has been checked using m.dot(m_inv).'
    print ''
    print '2) The matrix is not positive definite because one of its eigenvalues is a negative.'
    return m, m_inv,w,v, xlist

def get_condition_number(m):
    """
    http://mathworld.wolfram.com/ConditionNumber.html
    """
    u,s,v = np.linalg.svd(m);
    cond_num = max(s)/min(s)
    return cond_num

if __name__ == '__main__':
    #fname_pmbec = 'cov_matrix.pmbec.mat'
    #fname_pmbec_inv = 'InverseCovariance.txt' #This was used in the smmpmbec code.
    #m_inv = read_matrix(fname_pmbec_inv,header=False)
    #m = analyze_pmbec(offset=0.05)


    fname_blosum = 'ncbi.blosum62.sij.mat'
    result = analyze_blosum()
    m = result[0]
    m_inv = result[1]
    write_matrix(m, 'blosum62.mat')
    write_matrix(m_inv, 'blosum62.inv.mat')

    #Check that m_inv is really the inverse of 'm'.
    m_identity = m.dot(m_inv)
    #Set very small values to zeros.
    xsmall = 0.00001
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if abs(m_identity[i,j]) < xsmall:
                m_identity[i,j] = 0.0
            if abs(1.0 - m_identity[i,j]) < xsmall:
                m_identity[i,j] = 1.0
    plt.imshow(m_identity)
    plt.show()
    assert False not in (np.eye(20).flatten() == m_identity.flatten())
    




    
