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
#include "SMMSolve.h"

#include "BPException.h"
#include "math.h"

CSMMSolve::CSMMSolve(void)
{
}

CSMMSolve::~CSMMSolve(void)
{
}

void CSMMSolve::InitSolver(const CSMMSet &set, const InitParamSolve &p)
{
	m_tAA.SetToProduct(set.GetMatrix(),set.GetMatrix());
	m_y_inequal=set.GetMeasurements();	// m_y_inequal is used to store deviations from measurements for the inequality measurements.
	m_lambda_grouping = p.lambda_grouping;

	m_vec_length = m_tAA.NumCols();
	m_lambda_min = p.lambda_min;
	m_lambda_max = p.lambda_max;
	m_precision = p.precision;
	m_max_iterations = p.max_iterations;

	switch(m_lambda_grouping)
	{
	case ONE_LAMBDA:
		m_group_size=0;
		m_group_num=0;
		m_lambda.resize(1);
		m_inv_covar=CNumMat();
        break;
	case GROUP_LAMBDA:
		m_group_size = p.group_size;
		m_group_num  = (set.MatrixCols()-1)/m_group_size ;

		m_lambda.resize(m_group_num);
		m_inv_covar=CNumMat();
        break;
	case ONE_COVAR:
		m_group_size = p.group_size;
		m_group_num  = (set.MatrixCols()-1)/m_group_size ;
		m_lambda.resize(1);
		m_inv_covar=p.inverse_covar;
		assert(m_inv_covar.NumCols()==m_group_size);
        break;
	case GROUP_COVAR:
		m_group_size = p.group_size;
		m_group_num  = (set.MatrixCols()-1)/m_group_size ;
		m_lambda.resize(m_group_num);
		m_inv_covar=p.inverse_covar;
		assert(m_inv_covar.NumCols()==m_group_size);
		break;
	default:	throw BPException("Invalid lambda grouping");
	}
}


void CSMMSolve::SolveX(const CSMMSet &set, const CNumVec & param)
{
	assert(param.size()==m_lambda.size()); // assures that Init solver was called first

	// Convert param to lambda 
	const double param_max=log10(m_lambda_max);
	for(unsigned p=0; p<param.size(); p++)
	{
		if(param[p]>param_max)
			m_lambda[p]=m_lambda_min+m_lambda_max;
		else
			m_lambda[p]=m_lambda_min+pow(10.0,param[p]);
	}

	// Prepare tA x A + lambda (x Covariance)

	tAA_lam=m_tAA;
	switch(m_lambda_grouping)
	{
	case ONE_LAMBDA:
		for(unsigned int n=1; n<m_vec_length; n++)
			tAA_lam(n,n)+=m_lambda[0];
        break;
	case GROUP_LAMBDA:
		for(unsigned pos = 0; pos<m_group_num; pos++)
		{
			const unsigned offset = 1 + pos*m_group_size;
			for(unsigned i=0; i<m_group_size; i++)
				tAA_lam(offset+i,offset+i)+=m_lambda[pos];
		}
		break;
	case ONE_COVAR:
		for(unsigned pos = 0; pos<m_group_num; pos++)
		{
			const unsigned offset= 1 + pos * m_group_size;
			for(unsigned i=0; i<m_group_size; i++)
				for(unsigned j=0; j<m_group_size; j++)
					tAA_lam(offset+i,offset+j)+=m_inv_covar(i,j) * m_lambda[0];
		}
		break;
        break;
	case GROUP_COVAR:
		for(unsigned pos = 0; pos<m_group_num; pos++)
		{
			const unsigned offset= 1 + pos * m_group_size;
			for(unsigned i=0; i<m_group_size; i++)
				for(unsigned j=0; j<m_group_size; j++)
					tAA_lam(offset+i,offset+j)+=m_inv_covar(i,j) * m_lambda[pos];
		}
		break;
	default:	throw BPException("Invalid lambda grouping");
	}
	
	tAA_lam_inv.SetToInverse(tAA_lam);

	const CNumMat &A = set.GetMatrix();
	const CNumVec &y = set.GetMeasurements();
	CNumMat inverse;
	inverse.SetToProduct(tAA_lam_inv,A,false,true);

	if(set.InequalitiesPresent())
		CalcX_inequal(set, inverse);
	else
		m_x.SetToProduct(inverse,y);
}


void CSMMSolve::CalcX_inequal(const CSMMSet &set, const CNumMat & inverse)
{
	const CNumMat			&A		=set.GetMatrix();
	const CNumVec			&y		=set.GetMeasurements();
	const Vec<INEQUALITY>	&ineq	=set.GetInequalities();

	// Inequalities are present

	const unsigned	MAX_ITERATIONS	=	10000;
	const double	MAX_DIFFERENCE	=	1e-6;
	const double    ABSOLUTE_precision = m_precision*1e-3;

	m_x.SetToProduct(inverse,m_y_inequal);
	unsigned iteration=0;
	double relative_difference;
	
	CNumVec last_x;
	CNumVec ypred;

	do
	{
		ypred.SetToProduct(A,m_x);
		m_y_inequal=y;
		// Adjust m_y_inequal to ypred where allowed by inequality
		for(unsigned m=0; m<y.size(); m++)
		{
			if(ineq[m]!=EQUAL)
			{
				if(ineq[m]==GREATER)
				{
					if(ypred[m]>y[m])
						m_y_inequal[m]=ypred[m];
				}
				else if(ypred[m]<y[m])
					m_y_inequal[m]=ypred[m];
			}
		}

		last_x=m_x;
		m_x.SetToProduct(inverse,m_y_inequal);

		// Evaluate difference between new m_x and last_x
		relative_difference=0;
		for(unsigned n=0; n<m_x.size(); n++)
		{
			double max_abs=max(fabs(m_x[n]),fabs(last_x[n]));
			if(max_abs<ABSOLUTE_precision)
				continue;
			double rdiff=fabs(m_x[n]-last_x[n])/max_abs;
			if(rdiff>relative_difference)
				relative_difference=rdiff;
		}
		if(++iteration==MAX_ITERATIONS && relative_difference>MAX_DIFFERENCE)
		{
			BPNoConvergence e("No norm conversion for x=inverse(y_opt) with inequalities");
			e.m_message  << endl << "Iteration"		<< "\t" << MAX_ITERATIONS;
			e.m_message  << endl << "max difference"<< "\t" << relative_difference;
			e.m_message  << endl << "last x"		<< "\t"	<< last_x;
			e.m_message  << endl << "x"				<< "\t"	<< m_x;
			last_x-=m_x;
			e.m_message  << endl << "diff"			<< "\t"	<< last_x;
			e.m_message  << endl << "y"				<< "\t"	<< y;
			e.m_message  << endl << "ypred"			<< "\t"	<< ypred << endl;
			throw(e);
		}
	}
	while(relative_difference>MAX_DIFFERENCE);

	if(clog_detail.back()>DETAILED) 
		clog << endl << "Solve(lambda)\tIterations:\t" << iteration <<"\tm_x:\t" << m_x;
}


