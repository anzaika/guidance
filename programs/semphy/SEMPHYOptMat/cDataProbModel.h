// cDataProbModel.h: interface for the cDataProbModel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CDATAPROBMODEL_H__11E11909_3F96_45ED_9305_7335D9B04408__INCLUDED_)
#define AFX_CDATAPROBMODEL_H__11E11909_3F96_45ED_9305_7335D9B04408__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cProbModel.h"

class cDataProbModel {

public:
  inline const cOccuranceData& OccuranceData(void) const {return _D;};
	cDataProbModel(	const cOccuranceData& OccuranceData, 
					cProbModel& ProbModel);
	virtual ~cDataProbModel();

	cProbModel& ProbModel(void);

	cProbModel const& ProbModel(void) const;

	double NonConstLL(cSquareMatrix const& U) const;
	double LL(cSquareMatrix const& U) const;
    
	double NonConstLL(void) const;	
        double LL(void) const;

	cSquareMatrix DerivAllUxy(void) const;

	cRow DerivAllAngles(void);

	cRow DerivAllEigen(void) const;

	size_t size(void) const;

	double DerivRotxy(size_t i, size_t j) const;
	double DerivEigen(size_t i) const;
private:
	double MuabWhicht(const size_t a, const size_t b, const size_t Which_t) const;
  //	double DerivRotxy(size_t i, size_t j) const;
	double DerivUxy(size_t i, size_t j) const;
	cProbModel& _prob_model;
	const cOccuranceData& _D;
};

#endif // !defined(AFX_CDATAPROBMODEL_H__11E11909_3F96_45ED_9305_7335D9B04408__INCLUDED_)
