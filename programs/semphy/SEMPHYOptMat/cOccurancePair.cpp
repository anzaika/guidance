// cOccurancePair.cpp: implementation of the cOccurancePair class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cOccurancePair.h"
#include <string>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cOccurancePair::cOccurancePair(double Time, const cSquareMatrix& PseudoCounts) : 
	_time(Time), _pseudo_counts(PseudoCounts)
{
	if ( _time <= 0. ) {
		throw("cOccurancePair::cOccurancePair : non-positive time");
	}
}


void cOccurancePair::Read(istream& in)
{
  string tmp;
  in >> tmp;
  sscanf(tmp.c_str(),"dist=%lf",&_time);
  if ( _time <= 0. ) {
    throw("cOccurancePair::cOccurancePair : non-positive time");
  }
  _pseudo_counts.Read(in);
}


cOccurancePair::cOccurancePair(istream& in, const int length): _time(0.0), _pseudo_counts(length)
{
  this->Read(in);
}

// copy constructor
cOccurancePair::cOccurancePair(const cOccurancePair &Other) :
	_time(Other.Time()), _pseudo_counts(Other.Matrix())
{

}

cOccurancePair::~cOccurancePair(void)
{

}

cSquareMatrix cOccurancePair::NaiveNormalizeTime(void) const
{
	cSquareMatrix Result(_pseudo_counts);
	size_t i;
	for ( i = 0 ; i < Result.size() ; ++i ) {
		Result.Set(i,i,Result[i][i]-1.);
	}
	Result.ScaleSelf(1./_time);
	for ( i = 0 ; i < Result.size() ; ++i ) {
		Result.Set(i,i,Result[i][i]+1.);
	}
	return Result;
}

size_t cOccurancePair::size() const
{
	return _pseudo_counts.size();
}


double cOccurancePair::Time() const
{
	return _time;
}

const cSquareMatrix& cOccurancePair::Matrix() const
{
	return _pseudo_counts;
}

cOccurancePair cOccurancePair::operator =(const cOccurancePair &Other)
{
	_time = Other.Time();
	_pseudo_counts = Other.Matrix();
	return (*this);
}
 

void cOccurancePair::print(ostream& sout) const
{
  sout << "time="<< Time() << endl;
  _pseudo_counts.print(sout);
}
