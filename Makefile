COMPILER = g++
FLAGS = -lgsl -lgslcblas -lm
EXECUTABLE = smmpmbec

OBJECTS = BPException.o BPMath.o Dataset.o Letter.o MinimizerBase.o NumVec.o SMMCrossValidate.o SMMSeqPair.o SMMSet.o SMMSolve.o SeqMatrix.o SeqPair.o SeqSet.o VSet.o stdafx.o smm.o


$(EXECUTABLE) :  $(OBJECTS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) $(OBJECTS)

%.o : %.cpp
	$(COMPILER) -o $*.o -c $*.cpp

