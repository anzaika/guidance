// cCoeffMatrices.cpp: implementation of the cCoeffMatrices class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cCoeffMatrices.h"
#include "cDRotationMatrix.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cCoeffMatrices::cCoeffMatrices(const cAngles& Angles, int k)
  :	_k(k), _mleft(Angles.size()), _mright(Angles.size()), _angles(Angles)
{
  _mleft.Unit();
  _mright.Unit();
  // Nir fixed...
  _angles.RotateMatrix(_mright,_angles.NDims()-1,_k);
  _angles.RotateMatrix(_mleft,_k-1,-1);
}

cCoeffMatrices::~cCoeffMatrices()
{

}

cSquareMatrix cCoeffMatrices::CoeffMatrix(const cDRotationMatrix::tDrotType DrotType) const
{
	cAngles::tPlane CurrentPlane(_angles.GetDim(_k)); 

	cDRotationMatrix DRot(_angles.size(), CurrentPlane.first,CurrentPlane.second, DrotType);

	//	cout << "DROT"<<endl;
	//	DRot.print(cout);

	//	(_mleft*DRot.Matrix()*_mright).print(cout);

	return _mleft*DRot.Matrix()*_mright;
}
