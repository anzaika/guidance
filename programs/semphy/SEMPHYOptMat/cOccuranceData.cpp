// cOccuranceData.cpp: implementation of the cOccuranceData class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cOccuranceData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cOccuranceData::cOccuranceData(const int size): _mat_size(size)
{}

cOccuranceData::cOccuranceData(istream& in)
{
  int length, number;
  in >> length;
  in >> number;
  _mat_size=length;

  for (int i=0;i<number;++i){
    cOccurancePair p(in, length);
    AddPair(p);
  }
  _updateSuffisiontStatistics();
}


void cOccuranceData::_updateSuffisiontStatistics(void)
{
  _SstarDelta.resize(_mat_size);
  for (int a=0;a<_mat_size;++a){
    double s_a=0.0, sa_=0.0;
    for (int i=0;i<_data.size();++i) {
      for (int b=0;b<_mat_size;++b){
	sa_ += _data[i].Matrix()[a][b];
	s_a += _data[i].Matrix()[b][a];
      }
    }
    _SstarDelta[a]=s_a-sa_;
  }
}


cOccuranceData::cOccuranceData(const cOccurancePair& FirstPair)
{
	_data.push_back(FirstPair);
	_mat_size = FirstPair.size();		
  _updateSuffisiontStatistics();
}

cOccuranceData::~cOccuranceData()
{

}

cSquareMatrix cOccuranceData::AverageNormalization()
{
	cSquareMatrix SumAll(_mat_size);
	for ( tOccuranceData::iterator i = _data.begin() ; i != _data.end() ; ++i ){ 
		SumAll.AddSelf(i->NaiveNormalizeTime());
	}
	SumAll.ScaleSelf(1./_data.size()); // not _mat_size !!
	return SumAll;
}

void cOccuranceData::AddPair(const cOccurancePair &NextPair)
{
	_data.push_back(NextPair);
	if ( _mat_size != NextPair.size() ) {
		throw("cOccuranceData::AddPair : data of inappropriate matrix size");
	}
  _updateSuffisiontStatistics();
}

size_t cOccuranceData::alphabetSize() const
{
	return _mat_size;
}

// return S^(timeid)_[a,b] + S^(timeid)_[b,a]
double cOccuranceData::getSS(const int time_id,const int a,const int b) const 
{
  const cSquareMatrix& SS=Matrix(time_id);
  return (SS[a][b]+SS[b][a]);
}

const double& cOccuranceData::getSstar(const int a) const
{
  return (_SstarDelta[a]);
}

void cOccuranceData::print(ostream& sout) const
{
  int number=_data.size();
  sout<< number <<endl;
  for (tOccuranceData::const_iterator i=_data.begin();i!=_data.end();++i) {
    sout<<endl;
    i->print(sout);
  }
}
