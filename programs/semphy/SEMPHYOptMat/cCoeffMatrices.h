// cCoeffMatrices.h: interface for the cCoeffMatrices class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CCOEFFMATRICES_H__CD48CEE9_2163_4F9C_A6BA_8904DD86A547__INCLUDED_)
#define AFX_CCOEFFMATRICES_H__CD48CEE9_2163_4F9C_A6BA_8904DD86A547__INCLUDED_

#include "cOrthonormalMatrix.h"	// Added by ClassView
#include "cSquareMatrix.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cAngles.h"
#include "cDRotationMatrix.h"

class cCoeffMatrices  
{
public:
	int _k;
  cSquareMatrix CoeffMatrix(const cDRotationMatrix::tDrotType DrotType) const;
	cCoeffMatrices(const cAngles& Angles, int k);
	virtual ~cCoeffMatrices();

private:
	cOrthonormalMatrix _mleft;
	cOrthonormalMatrix _mright;
	const cAngles& _angles;

};

#endif // !defined(AFX_CCOEFFMATRICES_H__CD48CEE9_2163_4F9C_A6BA_8904DD86A547__INCLUDED_)
