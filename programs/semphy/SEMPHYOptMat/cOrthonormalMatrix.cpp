// cOrthonormalMatrix.cpp: implementation of the cOrthonormaMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cOrthonormalMatrix.h"
#include "cAngles.h"
#include <cmath>

using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cOrthonormalMatrix::cOrthonormalMatrix(size_t Size)
  : cSquareMatrix(Size,0.0)
{
	_initialize_diagonal();
}

cOrthonormalMatrix::~cOrthonormalMatrix(void)
{

}

void cOrthonormalMatrix::Set(size_t i, size_t j, double value)
{
	//Error("cOrthonormalMatrix::Set : method shaddowed");
}

void cOrthonormalMatrix::SetRow(size_t i, const cRow& row)
{
	//Error("cOrthonormalMatrix::SetRow : method shaddowed");
}

void cOrthonormalMatrix::Rotate(size_t i, size_t j, double sinus_alpha)
{
  double sin2a = sinus_alpha*sinus_alpha;
  if ( sin2a > 1 ) {
    throw("cOrthonormalMatrix::Rotate :  illegal sinus\n");
  }
  double cosinus_alpha = sqrt(1-sin2a);
  cRow NewRowI = GetRow(i)*cosinus_alpha;
  NewRowI.AddSelf(GetRow(j)*sinus_alpha);
  cRow NewRowJ = GetRow(j)*cosinus_alpha;
  NewRowJ.AddSelf(GetRow(i)*(-sinus_alpha));
  cSquareMatrix::SetRow(i,NewRowI);
  cSquareMatrix::SetRow(j,NewRowJ);
}

void cOrthonormalMatrix::RotateAlpha(size_t i, size_t j, double alpha)
{
  double sinus_alpha = sin(alpha);
  double cosinus_alpha = cos(alpha);
  cRow NewRowI = GetRow(i)*cosinus_alpha;
  NewRowI.AddSelf(GetRow(j)*sinus_alpha);
  cRow NewRowJ = GetRow(j)*cosinus_alpha;
  NewRowJ.AddSelf(GetRow(i)*(-sinus_alpha));
  cSquareMatrix::SetRow(i,NewRowI);
  cSquareMatrix::SetRow(j,NewRowJ);
}

cOrthonormalMatrix::cOrthonormalMatrix(const cProbs &BackgroundProbs)
  : cSquareMatrix(BackgroundProbs.size(),0.0)
{
  Unit();
  cAngles Angles(BackgroundProbs);
  Angles.RotateMatrix(*this);
}


// void cOrthonormalMatrix::print(ostream& sout)
// {
//   for (int i=0;i<_size;++i){
//     for (int j=0;j<_size;++j) {
//       sout << Get(i,j) <<" ";
//     }
//     sout <<endl;
//   }
// }
