#include "cRotationMatrix.h"
#include <cmath>

// eq:rotmat
cRotationMatrix::cRotationMatrix(const size_t size, const int i, const int j, const double phi)
  : cOrthonormalMatrix(size), _data(size)
{
  // the matrix is, at this stage,  eye(size)
	  _data.RotateAlpha(i,j,phi);
}

const cOrthonormalMatrix& cRotationMatrix::Matrix() const
{
	return _data;
}
