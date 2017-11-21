// cLineOpt.h: interface for the cLineOpt class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLINEOPT_H__851E9FAC_F197_45A4_959D_F3EC9120B533__INCLUDED_)
#define AFX_CLINEOPT_H__851E9FAC_F197_45A4_959D_F3EC9120B533__INCLUDED_

#include "cProbModelOptimizer.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cEigenVals.h"
#include "cOccuranceData.h"
#include "cProbs.h"
#include "cOrthonormalMatrix.h"

class cLineOpt  
{
public:
	void Opitmize(void);
	cLineOpt(const cOccuranceData& OccuranceData,
				   const cProbs& BackgroundProbs);
	virtual ~cLineOpt();
	const cEigenVals& Getd(void) const;
	const cOrthonormalMatrix& GetU(void) const;

private:
	double ImproveRot(double Delta);
	typedef void (*tOptFunc)(double);
	double ImproveEig(double Delta);
	cProbModelOptimizer _prob_model;
};


#endif // !defined(AFX_CLINEOPT_H__851E9FAC_F197_45A4_959D_F3EC9120B533__INCLUDED_)
