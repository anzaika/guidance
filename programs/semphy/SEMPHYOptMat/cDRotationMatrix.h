#ifndef __CDROTATIONMATRIX
#define __CDROTATIONMATRIX

#include "cSquareMatrix.h"
// The coefficient matrix Drot is not orthonormal, but just a matrix.
// it does not publicly inherit from cSquareMatrix in order to prevent non-const methods

  
class cDRotationMatrix 
{
public:
typedef enum { COS, SIN, ONE, NROTT } tDrotType;
	const cSquareMatrix& Matrix(void) const;


  cDRotationMatrix(const size_t size, const size_t i, const size_t j, const tDrotType t);
  virtual ~cDRotationMatrix(void);
  void print(ostream& sout = cout) const {_data.print(sout);};

private:
	cSquareMatrix _data;
};

#endif // __CDROTATIONMATRIX
