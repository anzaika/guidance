// cProbModel.h: interface for the cProbModel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CPROBMODEL_H__2DE37EF0_DAEA_4CE3_9897_7891B71392A0__INCLUDED_)
#define AFX_CPROBMODEL_H__2DE37EF0_DAEA_4CE3_9897_7891B71392A0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#ifdef _U
#undef _U
#endif
#ifdef _S
#undef _S
#endif

#include "cEigenVals.h"
#include "cOrthonormalMatrix.h"
#include "cProbs.h"
#include "cOccuranceData.h"

class cProbModelOptimizer  
{
public:
	const cEigenVals& Getd(void) const;
	const cOrthonormalMatrix& GetU(void) const;
	void AdjustAllEigen(double Delta);
	void RotateAll(double Delta);
	double NonConstLL(void) const;
	cProbModelOptimizer( const cOccuranceData& OccuranceData,
		const cProbs& BackgroundProbs);
	virtual ~cProbModelOptimizer(void);

private:
	double AvgDerivAllEigen(void) const;
	double DerivEigen(size_t i) const;
	double DerivRotxy(size_t i, size_t j) const;
	double DerivUxy(size_t i, size_t j) const;
	double Nabt(size_t i, size_t j, size_t Which_t) const;
	cEigenVals _d;
	cOrthonormalMatrix _U;
	const cOccuranceData& _S;
};

#endif // !defined(AFX_CPROBMODEL_H__2DE37EF0_DAEA_4CE3_9897_7891B71392A0__INCLUDED_)
