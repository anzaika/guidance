// cProbs.cpp: implementation of the cProbs class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cProbs.h"
#include <cmath>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cProbs::cProbs(size_t Size) : cRow(Size)
{
	(*this)[0] = 1.;
}

cProbs::~cProbs(void)
{

}

cProbs::cProbs(const cRow &row) : cRow(row)
{
  
	if ( fabs(Sum()-1.0) > CROW_EPSILON ){
		throw("cProbs::cProbs : vector does not sum up to 1");
	}
	if ( row < 0 ){
		throw("cProbs::cProbs : vector has negative elements");
	}
}

cProbs::cProbs(istream& in): cRow(in)
{
  if ( fabs(Sum()-1.0) > CROW_EPSILON ){
    throw("cProbs::cProbs : vector does not sum up to 1");
  }
  if ( (*this) < 0 ){
    throw("cProbs::cProbs : vector has negative elements");
  }
}
