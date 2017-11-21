// cRow.h: interface for the cRow class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CROW_H__D9B145FF_8F28_4C8D_BE26_6C6C5BA1792F__INCLUDED_)
#define AFX_CROW_H__D9B145FF_8F28_4C8D_BE26_6C6C5BA1792F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define CROW_EPSILON .0000001
#include <iostream>
#include<vector>
using namespace std;
typedef vector<double> tRow;
class cRow : public tRow  
{
public:
	void print(ostream& sout = cout) const;
	double Sum(void);
	cRow Exp(double Exponent) const;
	void ScaleSelf(double Scalar);
	void AddSelf(const cRow& Other);
	cRow operator+(const cRow& Other) const;
	cRow cRow::operator *(double Scalar) const;
	inline cRow& operator *=(const double Scalar);
	double Dot(const cRow& Other) const;
	cRow(size_t Size = 0);
	cRow(size_t Size, const double val);
	cRow(istream& in);
	virtual ~cRow(void);
	bool operator>=(double Scalar) const;
	bool operator<=(double Scalar) const;
	bool operator>(double Scalar) const;
	bool operator<(double Scalar) const;
	bool operator==(double Scalar) const;
	bool operator!=(double Scalar) const;
};

#endif // !defined(AFX_CROW_H__D9B145FF_8F28_4C8D_BE26_6C6C5BA1792F__INCLUDED_)
