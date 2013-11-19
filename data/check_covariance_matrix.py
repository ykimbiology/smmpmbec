

import numpy as np

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

def flatten(m):
    x=[]
    (nrow,ncol) = m.shape
    for i in range(nrow):
        for j in range(ncol):
            x.append(m[i,j])
    return x

if __name__ == '__main__':
    fname_pmbec = 'cov_matrix.pmbec.mat'
    fname_blosum = 'ncbi.blosum62.sij.mat'
    fname_pmbec_inv = 'InverseCovariance.txt'


    m = read_matrix(fname_pmbec,header=True) 
    m_inv_orig = read_matrix(fname_pmbec_inv,header=False) 

    m_inv = np.linalg.inv(m)

    mb = m + 0.04
    mb_inv = np.linalg.inv(mb)

    mc = m + np.eye(20)*0.02
    mc_inv = np.linalg.inv(mc)

    xa = flatten(m_inv_orig)
    xb = flatten(m_inv)
    xc = flatten(mb_inv)
    xd = flatten(mc_inv)

    plt.subplot(1,3,1)
    plt.plot(xa,xb, 'ko')

    # inv(M) == inv(M + c)?
    plt.subplot(1,3,2)
    plt.plot(xa,xc, 'ko')

    plt.subplot(1,3,3)
    plt.plot(xa,xd, 'ko')


    plt.show()


    
