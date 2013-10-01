#! /usr/bin/python

import sys
import os
import math
import cPickle
from xml.dom import minidom
import subprocess

class peptide:
    def __init__(self, seq="", meas = 0, ineq="=", pred = 0):
        assert(meas==float(meas))
        assert(ineq in "<=>")
        assert(pred==float(pred))
        for letter in seq:
            assert(letter in "ACDEFGHIKLMNPQRSTVWY")
        self.seq=seq
        self.meas=meas
        self.pred=pred
        self.ineq=ineq

class peptide_set:
    def __init__(self):
        self.set = []
        self.pep_length=0
        
    def load(self, in_file_name, pep_length=0):
        fin = open(in_file_name,"r")
        lines = fin.readlines()
        if pep_length<1:
            try:
                head = lines[0].split("\t")
                assert(len(head)==2)
                assert(int(head[1])==float(head[1]))
                self.pep_length= int(head[1])
                lines=lines[1:]
            except:
                print "Invalid header"
                print head
                raise
        else:
            self.pep_length
        for index, line in enumerate(lines):
            try:
                field=line.split()
                assert(len(field)==3)

                ineq = field[0].strip()
                meas = float(field[1])
                seq = field[2].strip()

                assert(len(seq)==self.pep_length)
                self.set.append(peptide(seq, meas,ineq))
            except:
                print "Exception"
                print index, field
                raise
        fin.close()

    def create_cv_sets(self, cv_num):
        sets = []
        for cv in range(cv_num):
            train = peptide_set()
            blind = peptide_set()
            train.pep_length=self.pep_length
            blind.pep_length=self.pep_length
            for index, pep in enumerate(self.set):
                if index % cv_num == cv:
                    blind.set.append(pep)
                else:
                    train.set.append(pep)
            sets.append([train, blind])
        return sets                                          

    def save_txt(self, out_name):
        out=open(out_name,"w")
        out.write("Seqlength:\t%d" %self.pep_length)
        for pep in self.set:
            out.write("\n%s\t%f\t%s\t%f" %(pep.ineq, pep.meas, pep.seq, pep.pred))
    def save_matlab(self, out_name):
        out=open(out_name,"w")
        for pep in self.set:
            out.write("1,")
            for letter in pep.seq:
                code = "ACDEFGHIKLMNPQRSTVWY".find(letter)
                assert(0<=code and code<20)
                for p in range(20):
                    if p == code:
                        out.write("1,")
                    else:
                        out.write("0,")
            out.write(pep.ineq +",")
            out.write(str(pep.meas) +"\n")
        out.close()
       
pep = peptide_set()
pep.load("RawData/murine-tap-9-mers_NO_EXTREME_OUTLIERS.txt")
pep.save_matlab("mathematica-peptides.csv")


