// cRow.cpp: implementation of the cRow class.
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "cRow.h"
#include <cmath>
using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cRow::cRow(size_t Size) : tRow(Size)
{
}

cRow::cRow(size_t Size, const double val) : tRow(Size, val)
{
}

cRow::~cRow()
{

}

void cRow::print (ostream& sout) const {
  for (int j=0;j<this->size();++j) {
    sout <<(*this)[j] <<" ";
  }
  sout <<endl;
}

cRow::cRow(istream& in){
  int length;
  in >> length;
  this->resize(length);
  for (int i=0;i<length;++i)
    in >> (*this)[i];
}



double cRow::Dot(const cRow &Other) const
{
  tRow::const_iterator i = begin();
  tRow::const_iterator j = Other.begin();
  double sum = 0;
  while( i != end() ) {
    sum += (*i)*(*j);
    ++i;
    ++j;
  }
  return sum;
}

cRow cRow::operator *(double Scalar) const
{
  cRow Result(*this);
  Result.ScaleSelf(Scalar);
  return Result;
}


inline cRow& cRow::operator *=(const double Scalar) 
{
  this->ScaleSelf(Scalar);
  return *this;
}

cRow cRow::operator +(const cRow &Other) const
{
  cRow Result(*this);
  Result.AddSelf(Other);
  return Result;
}

void cRow::AddSelf(const cRow &Other)
{
  tRow::iterator i = begin();
  tRow::const_iterator j = Other.begin();
  while( i != end() ) {
    (*i) += (*j);
    ++i;
    ++j;
  }
}

void cRow::ScaleSelf(double Scalar)
{
  for ( tRow::iterator i = begin() ; i != end() ; ++i ) {
    (*i) *= Scalar;
  }
}

cRow cRow::Exp(double Exponent) const
{
  cRow Result(size());
  tRow::const_iterator i = begin();
  tRow::iterator j = Result.begin();
  while( i != end() ) {
    (*j) = exp(*i*Exponent);
    ++i;
    ++j;
  }
  return Result;
}

double cRow::Sum(void)
{
  double sum = 0.;
  for ( tRow::iterator i = begin() ; i != end() ; ++i ) {
    sum += (*i);
  }
  return sum;
}
bool cRow::operator <=(double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) > Scalar) {
      return false;
    }
  }
  return true;
}

  bool cRow::operator >=(double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) < Scalar) {
      return false;
    }
  }
  return true;
}


bool cRow::operator ==(double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) != Scalar) {
      return false;
    }
  }
  return true;
}

  bool cRow::operator !=(double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) == Scalar) {
      return false;
    }
  }
  return true;
}
  


bool cRow::operator >(double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) <= Scalar) {
      return false;
    }
  }
  return true;
}

bool cRow::operator< (double Scalar) const
{
  for ( tRow::const_iterator i = begin() ; i != end() ; ++i ) {
    if ((*i) >= Scalar) {
      return false;
    }
  }
  return true;
}
  
