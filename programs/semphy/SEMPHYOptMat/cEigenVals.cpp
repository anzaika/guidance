// cEigenVec.cpp: implementation of the cEigenVec class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cEigenVals.h"

#include <cmath>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cEigenVals::cEigenVals(size_t Size) : _data(Size)
{
  _data[0]=0;
  for(int i=1;i<Size;++i)
    _data[i]=-1.0/(Size-1);
}

cEigenVals::~cEigenVals()
{
	
}

cEigenVals::cEigenVals(const cRow &row) : _data(row)
{
	_data[0] = 0;
	if ( fabs(_data.Sum() + 1)> CROW_EPSILON ) {
		// error
		throw("cEigenVals::EigenVals : illegal sum");
	}
}

void cEigenVals::IncreaseDecrease(size_t i, size_t j, double Offset)
{
	_data[i] += Offset;
	_data[j] -= Offset;
}

cRow cEigenVals::Exp(double t) const
{
	return _data.Exp(t);
}

const cRow& cEigenVals::Get(void) const
{
	return _data;
}


size_t cEigenVals::size(void) const
{
	return _data.size();
}
