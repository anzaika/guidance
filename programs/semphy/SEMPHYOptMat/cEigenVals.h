// cEigenVec.h: interface for the cEigenVec class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CEIGENVEC_H__08F19D3C_FAEC_411A_A215_91CA504F86EF__INCLUDED_)
#define AFX_CEIGENVEC_H__08F19D3C_FAEC_411A_A215_91CA504F86EF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cRow.h"
// maintain the invariants that the first element is 0, and the sum is -1
class cEigenVals
{
public:
	size_t size(void) const;
	inline double Get(size_t i) const;
        inline void Set(const int x, const double v )  { _data[x]=v;};
  cRow Exp(double t) const;
	const cRow& Get() const;
	void IncreaseDecrease(size_t i, size_t j, double Offset);
	cEigenVals(const cRow& row);
	cEigenVals(size_t Size);
	virtual ~cEigenVals();

private:
	cRow _data;
};

double cEigenVals::Get(size_t i) const
{
	return _data[i];
}

#endif // !defined(AFX_CEIGENVEC_H__08F19D3C_FAEC_411A_A215_91CA504F86EF__INCLUDED_)
