
This folder contains the two matrices of amino acid residue similarity: PMBEC and BLOSUM62.


1) cov_matrix.pmbec.mat
    Covariance matrix from the MHC Position Scanning Combinatorial Library binding data: PMBEC
    dimensions: 20x20

2) ncbi.blosum62.sij.mat
    NCBI's BLOSUM62 matrix. The matrix was used as if it was a covariance matrix in the bayesion prior model.
    dimensions: 20x20

3) InverseCovariance.txt
    Inverse of the PMBEC covariance matrix hardcoded in the c++ code.
    Note: A constant of 0.05 was added to all 20x20=400 elements of the matrix before taking the inverse.
    Note: This matrix is slightly different from the inverse of the matrix from 1). This difference appears to be due to slight differences in the PMBEC cov matrices.
    dimensions: 20x20