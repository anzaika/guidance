// cSquareMatrix.h: interface for the cSquareMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CSQUAREMATRIX_H__77D0D012_EE7E_403E_8FAA_1F9E39D07947__INCLUDED_)
#define AFX_CSQUAREMATRIX_H__77D0D012_EE7E_403E_8FAA_1F9E39D07947__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "cRow.h"
#include <iostream>
class cSquareMatrix  
{
public:
  inline const cRow& operator [](const int i) const {return _data[i];};
	double DotProduct(const cSquareMatrix& Other) const;
	void ScaleSelf(double Scalar);
	void AddSelf(const cSquareMatrix& Other);
	size_t size(void) const;
	cSquareMatrix Transpose(void) const;
	void TransposeSelf(void);
	inline void Set(size_t i, size_t j, double value) { _data[i][j] = value;};
	void SetRow(size_t i, const cRow& row);
	inline const cRow& GetRow(size_t i) const{return _data[i];};
	inline double Get(size_t i, size_t j) const {return _data[i][j];};
	typedef vector<cRow> tData;
	void Read(istream& in);
	cSquareMatrix  operator *(const cSquareMatrix &other) const;
	cSquareMatrix& operator *=(const cSquareMatrix &other);
	cSquareMatrix& operator *=(const double &Scalar) const { ScaleSelf(Scalar);return (this); }

	void print(ostream& sout = cout) const;
	cSquareMatrix(const size_t Size);
	cSquareMatrix(const size_t Size, const double val);
  // constract a diagonal matrix from the members of the diagonal
	cSquareMatrix(const cRow& diag);

	virtual ~cSquareMatrix(void);	
	void Unit(void);
	void _initialize_diagonal(void);
private:
	int _size;
	tData _data;
};

#endif // !defined(AFX_CSQUAREMATRIX_H__77D0D012_EE7E_403E_8FAA_1F9E39D07947__INCLUDED_)
