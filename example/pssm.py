#! /usr/bin/python

"""
Using the scoring matrix trained by smmpmbec, make binding predictions for a list of peptides.

INPUT: 
    1)Position Specific Scoring Matirx 
    2) A list of peptides
OUTPUT: 
    Predicted binding affinities for the peptides in log10(IC50).

"""

class PSSM:
    def __init__(self, fname):
        lines = open(fname,'r').readlines()
        lines_f = [l.strip() for l in lines if len(l.strip())>0]
        assert len(lines_f) == 22

        self.offset = float(lines_f[21].strip())
        self.mat={}
        for l in lines_f[1:21]:
            row = l.strip().split()
            res = row[0].strip()
            mat_row = map(float, row[1:])
            for pos in range(len(mat_row)):
                k = (pos, res)
                self.mat[k] = mat_row[pos]
        
    def predict_peptide(self, peptide):
        s = 0.0
        for pos, res in enumerate(peptide):
            k = (pos,res) # (int, character); e.g. (0, 'A')
            s = s + self.mat[k]
        s = s + self.offset
        return s

    def predict(self, peptide_list):
        score_list = [self.predict_peptide(p) for p in peptide_list]
        return score_list



if __name__ == '__main__':
    peptide_list = []
    peptide_list.append('AFAF')
    peptide_list.append('YSFL')
    peptide_list.append('AGAV')

    fname_pssm = 'output.txt'
    pssm = PSSM(fname_pssm)
    score_list = pssm.predict(peptide_list)
    header = ['peptide', 'predAffinity(log10(IC50))']
    print '\t'.join(header)
    for p,s in zip(peptide_list, score_list):
        print p, '\t', s
