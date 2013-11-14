
This folder contains the two matrices of amino acid residue similarity:

1) cov_matrix.pmbec.mat
    Covariance matrix from the mhc binding data: PMBEC
    dimensions: 20x20

2) ncbi.blosum62.sij.mat
    NCBI's BLOSUM62 matrix.
    dimensions: 20x20

3) InverseCovariance.txt
    Inverse of the PMBEC covariance matrix.
    Note: A constant of 0.05 was added to the matrix before taking the inverse.
    dimensions: 20x20