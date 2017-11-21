#ifndef __CROTATIONMATRIX
#define __CROTATIONMATRIX

#include "cOrthonormalMatrix.h"

class cRotationMatrix : public cOrthonormalMatrix
{
public:
	const cOrthonormalMatrix& Matrix(void) const;
  cRotationMatrix(const size_t size, const int i, const int j, const double phi);
  virtual ~cRotationMatrix(void){};
private:
	cOrthonormalMatrix _data;
	
};

#endif // __CROTATIONMATRIX
