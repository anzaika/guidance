// cSquareMatrix.cpp: implementation of the cSquareMatrix class.
//
//////////////////////////////////////////////////////////////////////
#include <vector>

#include "StdAfx.h"
#include "cSquareMatrix.h"

#include <cassert>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cSquareMatrix::cSquareMatrix(const size_t Size) : _size(Size),_data(0)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
		cRow Next(Size);
		_data.push_back(Next);
	}
}

cSquareMatrix::cSquareMatrix(const cRow& diag) : _size(diag.size()),_data(0)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
	  cRow Next(_size,0.0);
		Next[i]=diag[i];
		_data.push_back(Next);
	}
}

cSquareMatrix::cSquareMatrix(const size_t Size, const double val) : _size(Size),_data(0)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
	  cRow Next(Size, val);
		_data.push_back(Next);
	}
}

cSquareMatrix::~cSquareMatrix(void)
{

}

void cSquareMatrix::Read(istream& in)
{
  double val;
  for (int i=0;i<_size;++i)
    for (int j=0;j<_size;++j) {
      in >> val;
      Set(i,j,val);
    }
}

void cSquareMatrix::SetRow(size_t i, const cRow& row)
{
	_data[i] = row;
}

void cSquareMatrix::TransposeSelf(void)
{
	*this = Transpose();
}

void cSquareMatrix::Unit(void)
{
  ScaleSelf(0.0);
  _initialize_diagonal();
}


void cSquareMatrix::_initialize_diagonal(void)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
		this->Set(i,i,1.0);
	}
}


cSquareMatrix cSquareMatrix::Transpose(void) const
{
	cSquareMatrix Result(_size);
	for ( size_t i = 0 ; i < _size ; ++i ) {
		for ( size_t j = 0 ; j < _size ; ++j ) {
			Result._data[i][j] = this->_data[j][i];
		}
	}
	return Result;
}



size_t cSquareMatrix::size(void) const
{
	return _size;
}

void cSquareMatrix::AddSelf(const cSquareMatrix &Other)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
		(this->_data[i]).AddSelf(Other._data[i]);
	}
}

void cSquareMatrix::ScaleSelf(double Scalar)
{
	for ( size_t i = 0 ; i < _size ; ++i ) {
		(this->_data[i]).ScaleSelf(Scalar);
	}	
}

void cSquareMatrix::print(ostream& sout ) const
{
  for (int i=0;i<_size;++i)
    GetRow(i).print(sout);
  sout<<endl;
}


cSquareMatrix cSquareMatrix::operator *(const cSquareMatrix &other) const
{
  assert (other._size == _size);
  cSquareMatrix Result(this->_size,0.0);
  for (int i=0; i<_size; ++i)
    for (int j=0; j<_size; ++j)
      for (int k=0; k<_size; ++k)
	Result._data[i][j]+=this->_data[i][k]*other._data[k][j];
  return Result;
}


cSquareMatrix& cSquareMatrix::operator *=(const cSquareMatrix &other)
{
  assert (other._size == _size);
  cSquareMatrix Result(this->_size,0.0);
  for (int i=0; i<_size; ++i)
    for (int j=0; j<_size; ++j)
      for (int k=0; k<_size; ++k)
	Result._data[i][j]+=this->_data[i][k]*other._data[k][j];
  *this=Result;
  return(*this);
}


double cSquareMatrix::DotProduct(const cSquareMatrix &Other) const
{
  assert (Other._size == _size);
  double Result(0.0);
  for (int i=0; i<_size; ++i)
    for (int j=0; j<_size; ++j)
		Result+=this->_data[i][j]*Other._data[i][j];
  return(Result);
}
