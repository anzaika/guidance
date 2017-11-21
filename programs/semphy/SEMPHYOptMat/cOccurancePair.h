// cOccurancePair.h: interface for the cOccurancePair class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COCCURANCEPAIR_H__AEA486A9_E8DB_444F_8F2E_6C93C4C18D81__INCLUDED_)
#define AFX_COCCURANCEPAIR_H__AEA486A9_E8DB_444F_8F2E_6C93C4C18D81__INCLUDED_

#include "cSquareMatrix.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class cOccurancePair 
{
public:
	cOccurancePair(const cOccurancePair& Other);
	cOccurancePair(istream& in, const int length);
	const cSquareMatrix& Matrix(void) const;
	size_t size(void) const;
	// normalize by first order approximation.
	cSquareMatrix NaiveNormalizeTime(void) const;
	cOccurancePair(double Time, const cSquareMatrix& PseudoCounts);
	virtual ~cOccurancePair();
	double Time(void) const;
	cOccurancePair operator =(const cOccurancePair &Other);
	void print(ostream& sout = cout) const;
private:
	void Read(istream& in);
	double _time;
	cSquareMatrix _pseudo_counts;
};

#endif // !defined(AFX_COCCURANCEPAIR_H__AEA486A9_E8DB_444F_8F2E_6C93C4C18D81__INCLUDED_)
