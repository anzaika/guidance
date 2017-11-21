// cOrthonormalMatrix.h: interface for the cOrthonormalMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CORTHONORMAMATRIX_H__B22068D5_B34A_44DF_90EC_2FA8DFC11E74__INCLUDED_)
#define AFX_CORTHONORMAMATRIX_H__B22068D5_B34A_44DF_90EC_2FA8DFC11E74__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cSquareMatrix.h"
#include "cProbs.h"

class cOrthonormalMatrix : public cSquareMatrix  
{
public:
	cOrthonormalMatrix(const cProbs& BackgroundProbs);
	void RotateAlpha(size_t i, size_t j, double alpha);
	void Rotate(size_t i, size_t j, double sinus_alpha);
  //        void print(ostream& sout = cout);
	cOrthonormalMatrix(size_t size);
	virtual ~cOrthonormalMatrix(void);
	
protected:
	// hide set operations, which are now longer legal
	void Set(size_t i, size_t j, double value);
	void SetRow(size_t i, const cRow& row);
};

#endif // !defined(AFX_CORTHONORMAMATRIX_H__B22068D5_B34A_44DF_90EC_2FA8DFC11E74__INCLUDED_)
