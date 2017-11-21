#include "cDRotationMatrix.h"

cDRotationMatrix::cDRotationMatrix(const size_t size, const size_t i, 
								   const size_t j, const tDrotType t)
: _data(size,0)
{
  switch ( t ) {
  case ONE:
    for ( int k = 0 ; k < size ; ++k ) {
      _data.Set(k,k,1.);
    }
    _data.Set(i,i,0.);
    _data.Set(j,j,0.);
    break;
  case COS:
    _data.Set(i,i,1.);
    _data.Set(j,j,1.);
    break;
  case SIN:
    _data.Set(i,j,1.);
    _data.Set(j,i,-1.);
    break;
  default:
    throw("cDRotationMatrix::cDRotationMatrix : unknown type");
  }
}

cDRotationMatrix::~cDRotationMatrix()
{
}

const cSquareMatrix& cDRotationMatrix::Matrix() const
{
	return _data;
}

