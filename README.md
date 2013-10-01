

INTRODUCTION
============
The algorithm SMMPMBEC models binding interactions between MHC class I molecules and peptide ligands.
Given a list of (peptide, affinity), it returns a PSSM capturing binding energy contributions of residues at various positions.


REFERENCE
=========
Kim Y, Sidney J, Pinilla C, Sette A, and Peters B. "Derivation of an amino acid similarity matrix for peptide:MHC binding and its application as a Bayesian prior" BMC Bioinformatics 2007.
http://www.ncbi.nlm.nih.gov/pubmed/19948066


REQUIREMENTS
============
1) linux environment
2) GSL [http://www.gnu.org/software/gsl/]

COMPILING
=========
Simply issue the command 'make' within the folder.
This should generate an executable 'smmpmbec'.


TRAINING
========
See ./example/run_smmpmbec.sh for a training example.
