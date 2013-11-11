

Introduction
============
The algorithm SMMPMBEC models binding interactions between MHC class I molecules and peptide ligands.
Given a list of (peptide, affinity), it returns a PSSM capturing binding energy contributions of residues at various positions.


Reference
=========
Kim Y, Sidney J, Pinilla C, Sette A, and Peters B. "Derivation of an amino acid similarity matrix for peptide:MHC binding and its application as a Bayesian prior" BMC Bioinformatics 2007.
http://www.ncbi.nlm.nih.gov/pubmed/19948066


Requirements
============
1) Linux environment with a gcc compiler (e.g. Ubuntu).
2) Gnu Scientific Library(GSL) [http://www.gnu.org/software/gsl/]

Compiling
=========
Simply issue the command 'make' within the folder.
This should generate an executable 'smmpmbec'.


Training
========
See ./example/run_smmpmbec.sh for a training example.
Trained model is a position specific scoring matrix (PSSM).
For binding data set for peptides of 9 residues, the matrix has dimensions of 20 residues by 9 positions: 20x9.

Making predictions
==================
For a given peptide, its predicted binding affinity in log10(IC50) is a sum of corresponding entries in the PSSM.
