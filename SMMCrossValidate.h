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
#pragma once

#include "SMMSolve.h"
#include "MinimizerBase.h"

// Calculates the cross validated distance of predictions to measurements 
// for a given lambda 


struct InitParamCV : public InitParamSolve
{
	unsigned srand;
	unsigned train_repeats;
	unsigned tolerated_non_convergence;
	unsigned cv_num;
	double lambda_start;

	InitParamCV()
	{
 		srand = 0;
		train_repeats=10;
		tolerated_non_convergence = 3;
		cv_num =10;
		lambda_start = 1.0;
	};

	void log(ostream &out) const
	{
		out << endl << "CV Param:";
		out << endl << "srand\t" << srand;
		out << endl << "train_repeats\t" << train_repeats;
		out << endl << "tolerated_non_convergence\t" << tolerated_non_convergence;
		out << endl << "cv_num\t" << cv_num;
		out << endl << "lambda_start\t"	<< lambda_start;
		InitParamSolve::log(out);
	}
};


class CSMMCrossValidate:   private CMinimizerBase 
{
public:
	CSMMCrossValidate(void) {};
	~CSMMCrossValidate(void) {};

public:	
	CNumVec	Train(const CSMMSet &train_set, const InitParamCV & cvp);
protected:
	CNumVec	_Train(const CSMMSet &train, const InitParamCV & cvp);
	virtual double Distance(const CNumVec & param);
	virtual double Gradient(const CNumVec & param, CNumVec &p_gradient);
private:
	double			m_precision;
	Vec<CSMMSolve>	m_solver;
	CVSet<CSMMSet>  m_train_set;	
	CNumVec gradient_current; // YK: To keep record gradient of cv distance function.
};
