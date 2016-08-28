// cProbModel.cpp: implementation of the cProbModel class.
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "cProbModelOptimizer.h"

#include <cmath>

using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
	
cProbModelOptimizer::cProbModelOptimizer(const cOccuranceData& OccuranceData,
										 const cProbs& BackgroundProbs) :
	_d(BackgroundProbs.size()), _U(BackgroundProbs), _S(OccuranceData)
{
	if ( BackgroundProbs.size() != OccuranceData.alphabetSize() ) {
		throw("cProbModelOptimizer::cProbModelOptimizer : size mismatch");
	}
}

cProbModelOptimizer::~cProbModelOptimizer()
{

}

double cProbModelOptimizer::Nabt(size_t i, size_t j, size_t Which_t) const
{
	cRow Expt = _d.Exp(_S.Time(Which_t));
	double sum = 0.;
	for ( size_t k = 0 ; k < _U.size() ; ++k ) {
		sum += Expt[k]*_U[k][i]*_U[k][j];
	}
	return sum;
}

double cProbModelOptimizer::NonConstLL(void) const
{
	double sum = 0.;
	for ( size_t Which_t = 0 ; Which_t < _S.size() ; ++Which_t ) {
		const cSquareMatrix CurrentMat(_S.Matrix(Which_t));
		for ( size_t i = 0 ; i < _U.size() ; ++i ) {
			for ( size_t j = 0 ; j < _U.size() ; ++j ) {
				sum += CurrentMat[i][j]*log(Nabt(i,j,Which_t));
			}
		}
	}
	return sum;
}

double cProbModelOptimizer::DerivUxy(size_t i, size_t j) const
{
	double sum = 0.;
	//	double dx = _d.Get(i);
	for ( size_t Which_t = 0 ; Which_t < _S.size() ; ++Which_t ) {
	  //		double ExpVal = _S[Which_t].Time()*dx;
		//MN		double etdx = exp(ExpVal);
		const cSquareMatrix CurrentMat(_S.Matrix(Which_t));
		for ( size_t k = 0 ; k < _U.size() ; ++k ) {
			sum += _U[i][k]*(CurrentMat[k][j]+CurrentMat[j][k])
				/Nabt(i,j,Which_t);
		}
	}
	return sum;
}

double cProbModelOptimizer::DerivRotxy(size_t i, size_t j) const
{
	double sum = 0.;
	for ( size_t k = 0 ; k < _U.size() ; ++k ) {
		sum += DerivUxy(i,k);
		sum -= DerivUxy(j,k);
	}
	return sum;
}


double cProbModelOptimizer::DerivEigen(size_t i) const
{
	double sum = 0.;
	for ( size_t Which_t = 0 ; Which_t < _S.size() ; ++Which_t ) {
		const cSquareMatrix CurrentMat(_S.Matrix(Which_t));
		double t = _S.Time(Which_t);
		double etdx = exp(t*_d.Get(i));
		for ( size_t j = 0 ; j < _U.size() ; ++j ) {
			for ( size_t k = 0 ; k < _U.size() ; ++k ) {
				sum += CurrentMat[j][k]*t*etdx*_U[i][j]*_U[i][k]
					/Nabt(j,k,Which_t);
			}
		}
	}
	return sum;
}


double cProbModelOptimizer::AvgDerivAllEigen(void) const
{
	double sum = 0.;
	for ( size_t i = 1 /* first eigenvalue is constant */ ; i < _d.size() ; ++i ) {
		sum += DerivEigen(i);
	}
	return sum/(_d.size()-1);
}

void cProbModelOptimizer::RotateAll(double Delta)
{
	for ( size_t i = 1 ; i < _U.size()-1 ; ++i ) {
		for ( size_t j = i+1 ; j < _U.size() ; ++j ) {
			_U.Rotate(i,j,Delta*DerivRotxy(i,j));
		}
	}
}

void cProbModelOptimizer::AdjustAllEigen(double Delta)
{
	double AvgGrad = AvgDerivAllEigen();
	double NeedToAdd = 0;
	for ( size_t i = 1 ; i < _d.size()-1 ; ++i ) {
		double NewAdd = DerivEigen(i)-AvgGrad+NeedToAdd;
		_d.IncreaseDecrease(i,i+1, NewAdd);
		NeedToAdd = NewAdd;
	}
}

const cOrthonormalMatrix& cProbModelOptimizer::GetU(void) const
{
	return _U;
}

const cEigenVals& cProbModelOptimizer::Getd(void) const
{
	return _d;
}
