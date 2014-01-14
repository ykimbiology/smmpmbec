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
#include "SMMCrossValidate.h"
#include <algorithm>
#include "math.h"



CNumVec	CSMMCrossValidate::Train(const CSMMSet &train_set, const InitParamCV & init_param)
{
	// Repeats the _Train process train_repeats times, and averages the solutions.
	// In each iteration, the cross validation is done on different sets, which
	// makes the average solution more robust.

	init_param.log(clog);
	srand(init_param.srand);

	CNumVec mean_solution(train_set.MatrixCols());
	mean_solution = 0;
	unsigned count_solution=0;
	unsigned non_convergence=0;

	bool isConverged = true;

	for(unsigned r=0; r<init_param.train_repeats; r++)
	{
		try
		{
			mean_solution+= _Train(train_set, init_param);
			count_solution++;
		}
		catch(BPNoConvergence b)
		{
			non_convergence++;
			clog << endl << "Repeat:\t" << r  <<"\t" << non_convergence <<"\t failure to converge. Message: " << endl;
			clog << b.m_message.str();
			if(non_convergence > init_param.tolerated_non_convergence)
				isConverged = false;
				throw BPNoConvergence("Multiple failures to converge");
		}
	}

	if (isConverged == true) {

		mean_solution*=1.0/count_solution;
		if(clog_detail.back()>=MEDIUM)
		clog << endl << "Final:\t\t" << mean_solution;
		cout << "model_final\t" << mean_solution << endl;

		return(CNumVec(mean_solution));
	} else {
		mean_solution.clear();
		return(CNumVec(mean_solution));
	}
}


CNumVec CSMMCrossValidate::_Train(const CSMMSet &train, const InitParamCV & init_param)
{
	const unsigned & cv_num = init_param.cv_num;
	m_precision = init_param.precision;

	// Splits the training set into cv_num cross validation sets
	// Sets up solvers for each of these
	m_train_set = train;
	m_train_set.GenerateSubsets(cv_num);
	m_solver.resize(cv_num);

	// In all cases, start with ONE_LAMBDA or ONE_COVAR to determine a good start lambda
	InitParamCV init_param2 = init_param;
	switch(init_param.lambda_grouping)
	{
	case GROUP_LAMBDA: init_param2.lambda_grouping=ONE_LAMBDA;
	case GROUP_COVAR:  init_param2.lambda_grouping=ONE_COVAR;
	}

	// Set up solvers
	for(unsigned c=0; c<cv_num; c++)
		m_solver[c].InitSolver(m_train_set.SubsetTrain(c), init_param2);
	CNumVec param(1);
	param=log(init_param.lambda_start);
	// Optimize
	double dist = MinimizeSteepestDescentOneDimensional(param,m_precision);

	if(init_param.lambda_grouping!=init_param2.lambda_grouping)
	{
		// Now do complete training
		// Set up solvers
		for(unsigned c=0; c<cv_num; c++)
			m_solver[c].InitSolver(m_train_set.SubsetTrain(c), init_param);
		// Start with previously determined parameters
		double para0 = param[0];
		param.resize(m_solver[0].GetLambda().size());
		param=para0;
		// Optimize
		dist = MinimizeSteepestDescent(param,m_precision);
		cout << "cv_dist\t" << dist << endl;
	}
	// Calculate the solution as an average of the individual cross validation solutions
	CNumVec mean_solution(m_solver[0].GetX().size());
	mean_solution = 0;
	for (unsigned c=0; c<cv_num; c++) {
		mean_solution += m_solver[c].GetX();
		cout << "x_i\t" << m_solver[c].GetX() << endl;
	}

	CNumVec lambdas_log(m_solver[0].GetLambda().size());
	lambdas_log = 0;
	for (int i=0; i<lambdas_log.size(); i++) {
		lambdas_log[i] = log(m_solver[0].GetLambda()[i]);  // In log() space; log="natural log".
	}
	//cout << "lambdas_log\t" << m_solver[0].GetLambda() << endl; // YK: For one 10-fold cv, there is one set of lambdas optimized.
	cout << "lambdas_log\t" << lambdas_log << endl;

	cout << "dfdlambdas_log\t" << gradient_current << endl;

	mean_solution*=1.0/cv_num;
	cout << "x_mean\t" << mean_solution << endl;

	if(clog_detail.back()>=MEDIUM)
		clog << endl << "Dist:\t" << dist << "\tLambda:\t" << m_solver[0].GetLambda() << "Solution:\t" << mean_solution;
	return (CNumVec(mean_solution));
}


/**
 * INPUT: param, regularization parameters.
 * OUTPUT: Returns an internal cross-validated distance between predicted and measured affinities.
 *
 */
double CSMMCrossValidate::Distance(const CNumVec &param)
{
	double dist=0;
	unsigned count=0;
	for(unsigned c=0; c<m_train_set.SubsetNum(); c++)
	{
		const CSMMSet &train=m_train_set.SubsetTrain(c);
		CSMMSet		  &blind=m_train_set.SubsetBlind(c);
		CSMMSolve	  &solver=m_solver[c];

		solver.SolveX(train,param);
		dist +=blind.XDistance(solver.GetX());
		count+=blind.ElementNumber();
	}
	dist/=count; // YK: Mean of sum of squared differences between predicted and measured in log space.

	if(clog_detail.back()>DETAILED)
		clog << endl << "Dist\t"  << dist << "\tLambda\t"  << param;
	return dist;
}



double CSMMCrossValidate::Gradient(const CNumVec &param, CNumVec &gradient)
{
	gradient.resize(param.size());
	double dist0=Distance(param);
	double delta = m_precision / 1000.0;
	for (unsigned p=0; p<param.size(); p++)
	{
		CNumVec param2=param;
		param2[p]+=delta;
		double dist=Distance(param2);
		gradient[p]=(dist-dist0)/delta; // YK: Finite difference technique to calculate the gradient.
		gradient_current = gradient;    // YK: Keep a record of the current gradient.
	}
	return(dist0);


	// Could be made analytically, but is unlikely to be much quicker, and most likely harder to maintain.

	/*
	assert(gradient.size()==param.size());

	double dist=0;
	unsigned count=0;
	gradient=0;

	for(unsigned c=0; c<m_train_set.SubsetNum(); c++)
	{
		const CSMMSet &train=m_train_set.SubsetTrain(c);
		CSMMSet		  &blind=m_train_set.SubsetBlind(c);
		CSMMSolve	  &solver=m_solver[c];

		solver.SolveXdLambda(train,lambda,m_precision, m_max_norm_iterations);
		const CNumMat &dxdl=solver.Getdxdlambda();
		const CNumVec &x   =solver.GetX();

		CNumVec ddistdx;
		dist+=blind.XGradient(x,ddistdx);
		count+=blind.ElementNumber();

		CNumVec g;
		g.SetToProduct(dxdl,ddistdx);
		gradient+=g;
	}
	gradient*=1.0/count;
	dist/=count;
	if(clog_detail.back()>DETAILED)
	{
		clog << endl << "Dist\t"  << dist << "\tLambda\t"  << lambda;
		clog << endl << "Dist\t"  << dist << "\tGrad\t"  << gradient;
	}
	//	return dist;

	assert(false);
// Here the Gradient Conversion starts
// needs morefixing
//	CalcLambdaFromParameters(param);
//	CNumVec l_gradient(m_lambda.size());
//	double dist=CVGradient(l_gradient,m_lambda);
	CNumVec l_gradient = gradient;

	p_gradient.resize(param.size());
	for(unsigned p=0; p<param.size(); p++)
	{
		const UIntVec & lambda_group=m_lambda_groups[p];
		p_gradient[p]=0;
		for(unsigned i=0; i<lambda_group.size(); i++)
			p_gradient[p]+=m_lambda[lambda_group[i]]*l_gradient[lambda_group[i]];
	}

	if(m_only_lambda_descent && p_gradient.size()>1)	// Disallowing anything but descent for the lambda values actually makes the search for a minimum more robust - for a starting point found with FindStartLambda
	{
		for(unsigned n=0; n<p_gradient.size(); n++)
		{
			if(p_gradient[n]<0)
				p_gradient[n]=0;
		}
	}

	if(clog_detail.back()>MEDIUM)
	{
		clog << endl << "Dist\t"  << dist << "\tParam\t"  << param;
		clog << endl << "Dist\t"  << dist << "\tGrad\t"  << p_gradient;
	}
	return(dist);
	*/
}



/*
void   CSMMCrossValidate::AdjustParameters(CNumVec &param)
{
	const double param_max=log(m_lambda_max);
	for(unsigned p=0; p<param.size(); p++)
	{
		if(param[p]>param_max)
			param[p]=param_max;
	}
}

void	CSMMCrossValidate::CreateMaskedCopy(const CSMMCrossValidate & cv, const UIntVec & mask)
{
	m_train_set.CreateMaskedCopy(cv.m_train_set,mask);
	m_solver.resize(m_train_set.SubsetNum());
	for(unsigned c=0; c<m_train_set.SubsetNum(); c++)
		m_solver[c].InitSolver(m_train_set.SubsetTrain(c), m_use_covar_matrix, m_inv_covar_matrix);
}

*/
