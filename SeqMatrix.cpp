/*
www.mhc-pathway.net/smm
Original file by Bjoern Peters.

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any
damages arising from the use of this software.

Permission is granted to anyone to use this software for any
purpose, including commercial applications, and to alter it and
redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and
must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.
*/
#include "stdafx.h"
#include "SeqMatrix.h"
#include "BPMath.h"




void CSeqMatrix::SMMTrain(const CSeqSet &seqset, InitParamCV init_param)
{
	clog << endl << "CSeqMatrix::SMMTrain(" << seqset.GetName()<<")";
	clog << endl << "Elements in Set: " << seqset.ElementNumber();
	resize(seqset.SeqLength());
	CSMMSet smm_train;
	seqset.ConvertToSMMSet(smm_train,true);
	CSMMCrossValidate smm;
	init_param.group_size=CLetter::AlphabetSize();
	m_matrix = smm.Train(smm_train, init_param);
}


void CSeqMatrix::Load(const string & name)
{
	ex_ifstream in(name);
	try
	{
		int length;
		string x;
		in >> x >> length;
		if(in.fail())
			throw(BPException("Could not load NumCols for CSeqMatrix"));
		Load(in,length,false);
	}
	catch(BPException e)
	{
		throw BPFileException(e,name);
	}
}




void CSeqMatrix::Save(const string & name)	const
{
	ex_ofstream out(name);
	out <<"NumCols:\t" << NumCols() <<endl;
	for(CSequence::const_iterator let=CLetter::Alphabet().begin(); let<CLetter::Alphabet().end(); let++)
	{
		out << let->Letter();
		for(unsigned pos=0; pos<NumCols(); pos++)
			out << "\t" << Matrix(*let,pos);
		out <<endl;
	}
	out << endl << endl << Offset();
}

void CSeqMatrix::Predict(CSeqSet  &set)	const
{
	const Vec<CSequence> &seqs=set.GetSequences();
	CNumVec predictions(seqs.size());
	for(unsigned i=0; i<seqs.size(); i++)
		predictions[i]=Score(seqs[i]);
	set.Predict(predictions);
};


double CSeqMatrix::Score(CSequence::const_iterator it, CSequence::const_iterator end)	const
{
	assert(end-it==NumCols());
	double score=Offset();
	for(unsigned i=0; i<NumCols(); i++)
		score+=Matrix(*(it+i),i);
	return(score);
}

void CSeqMatrix::Load(const string & name, int length, bool no_offset)
{
	ex_ifstream in(name);
	try
	{
		Load(in,length,no_offset);
	}
	catch(BPException e)
	{
		throw BPFileException(e,name);
	}

}

void CSeqMatrix::Load(istream & in, int length, bool no_offset)
{
	resize(length);
	unsigned num=0;
	while(!in.fail() && num<CLetter::AlphabetSize())
	{
		CLetter let;
		in >> let;
		num++;
		for(int pos=0; pos<length; pos++)
			in >> Matrix(let,pos);
		if(in.fail())
			throw BPException("Invalid CSeqMatrix file");
	}
	if(num!=CLetter::AlphabetSize())
		throw BPException("Incomplete CSeqMatrix file");
	if(no_offset)
		Offset()=0.0;
	else
	{
		in >> Offset();
		if(in.fail())
			throw BPException("No offset in CSeqMatrix file");
	}
}

unsigned CSeqMatrix::GetMatrixSize() {

	//cout << "CSeqMatrix::IsMatrixEmpty " << m_matrix.size() << endl;
	return m_matrix.size();
}

void CSeqMatrix::AdjustMatrix(CSeqSet & set)
{
	CNumVec meas=set.GetMeasurements();
	Vec<INEQUALITY> ineq=set.GetInequalities();

	Predict(set);
	double m=1,b=0;
	double diff_norm;
	CNumVec last_pred=set.GetPredictions();
	do
	{
		LinearRegression(set.GetPredictions(),meas,m,b);
		CNumVec pred=set.GetPredictions();
		pred*=m;
		pred+=b;
		for(unsigned n=0; n<pred.size(); n++)
		{
			if(pred[n]>meas[n] &&ineq[n]==GREATER)
				meas[n]=pred[n];
			else if(pred[n]<meas[n] &&ineq[n]==SMALLER)
				meas[n]=pred[n];
		}
		CNumVec diff=last_pred;
		diff-=pred;
		last_pred=pred;
		diff_norm=diff.ScalarProduct(diff);
	}
	while(diff_norm>.000001);
	m_matrix*=m;
	m_matrix[0]+=b;
}


void CSeqMatrix::AdjustOffset(CSeqSet & set)
{
	CNumVec meas=set.GetMeasurements();
	Vec<INEQUALITY> ineq=set.GetInequalities();

	Predict(set);
	double diff_norm;
	CNumVec last_pred=set.GetPredictions();

	double mean_diff=0;
	do
	{
		CNumVec pred=set.GetPredictions();
		mean_diff=0;
		for(unsigned n=0; n<pred.size(); n++)
			mean_diff+=meas[n]-pred[n];
		mean_diff/=pred.size();
		pred+=mean_diff;
		for(unsigned n=0; n<pred.size(); n++)
		{
			if(pred[n]>meas[n] &&ineq[n]==GREATER)
				meas[n]=pred[n];
			else if(pred[n]<meas[n] &&ineq[n]==SMALLER)
				meas[n]=pred[n];
		}
		CNumVec diff=last_pred;
		diff-=pred;
		last_pred=pred;
		diff_norm=diff.ScalarProduct(diff);
	}
	while(diff_norm>.000001);
	m_matrix[0]+=mean_diff;
}




void CSeqMatrix::AverageMatrix(const Vec<CSeqMatrix> matrices)
{
	resize(matrices.front().NumCols());
	assert(matrices.size()>0);
	m_matrix=matrices.front().m_matrix;
	for(unsigned n=1; n<matrices.size(); n++)
		m_matrix+=matrices[n].m_matrix;
	m_matrix*=1.0/matrices.size();
}

void CSeqMatrix::CountSet(const CSeqSet & set)
{
	resize(set.SeqLength());
	const Vec<CSequence> &seqs=set.GetSequences();
	for(Vec<CSequence>::const_iterator s=seqs.begin(); s<seqs.end(); s++)
		for(unsigned pos=0; pos<NumCols(); pos++)
			Matrix((*s)[pos],pos)++;
}


