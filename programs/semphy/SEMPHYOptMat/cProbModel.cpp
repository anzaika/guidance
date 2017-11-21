//////////////////////////////////////////////////////////////////////
// cProbModel.cpp: implementation of the ProbModel class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cProbModel.h"
#include "cOptManager.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


cProbModel::cProbModel(const cProbs& BackgroundProbs) :
	_d(BackgroundProbs.size()), _U(BackgroundProbs.size()), _Phi(BackgroundProbs)
{
  _updateU();
}

cProbModel::~cProbModel()
{

}

double cProbModel::Muabt(const size_t a, const size_t b, const double t) const
{

//   cout <<"t="<<t<<"   d:";
//   _d.Get().print(cout);
//   cout <<endl<<"EXP(t*d)=";
  cRow Expt = _d.Exp(t);
//   Expt.print(cout);
//   cout <<endl;
  double sum = 0.;
  for ( size_t k = 0 ; k < _U.size() ; ++k ) {
    sum += _U[k][a]*Expt[k]*_U[k][b];
  }
  return sum;
}



cSquareMatrix cProbModel::Mut(const double t) const
{
  cSquareMatrix Diag(_d. Exp(t));	// diagonal matrix with exp(t*d)

//   _U.print();
//   Diag.print();
  
  return _U.Transpose()*Diag*_U;
}

// give a diagonal matrix and get the resulting mut
cSquareMatrix cProbModel::MutD(const double t, const cRow& D) const
{
  cSquareMatrix Diag(D. Exp(t));	// diagonal matrix with exp(t*D)

  //  cout<< "MutD"<<endl;
  // _U.print(cout);
  // Diag.print(cout);
  
  return _U.Transpose()*Diag*_U;
}


const cOrthonormalMatrix& cProbModel::GetU(void) const
{
	return _U;
}

const cEigenVals& cProbModel::GetD(void) const
{
	return _d;
}
const cAngles& cProbModel::Angles() const
{
	return _Phi;
}

void cProbModel::setAngles(const cAngles& newAngles)
{
  //  cout << "SetAngles"<< endl;
  for (int k=_Phi.size()-1;k<_Phi.NDims();++k){

    //        cout << _Phi.GetDim(k).first<<"&"<<_Phi.GetDim(k).second<<endl;
    _Phi.Set(k,newAngles.get(k));
  }
  //    _Phi.print(cout);
  _updateU();
}

void cProbModel::_updateU(void)
{
  _U.Unit();
  _Phi.RotateMatrix(_U);
}


void cProbModel::setAngle(const int k, const double newAngle)
{
  _Phi.Set(k,newAngle);
  _updateU();
}


// by default, t=1.0
cSquareMatrix cProbModel::GetQ() const
{
  cSquareMatrix Diag(_d.Get());	// diagonal matrix with exp(t*d)

  return _U.Transpose()*Diag*_U;
}

cSquareMatrix cProbModel::GetR() const
{
  cSquareMatrix psq(GetPiSqrt());	
  cSquareMatrix psqr(GetPiSqrtRev());	

  cout<< "GetPiSqrt";
  GetPiSqrt().print();
  cout <<endl;
  cout <<endl;
  psq.print();
  cout <<endl;
  cout<< "GetPiSqrtRev";
  GetPiSqrtRev().print();
  cout <<endl;
  cout <<endl;
  psqr.print();
  cout <<endl;
  
  cout<< "Q"<<endl;
  GetQ().print();
  cout <<endl;
  
  cout<< "R shoud be"<<endl;
  (psqr*GetQ()*psq).print();

  return psqr*GetQ()*psq;
}
  
cProbs cProbModel::GetPi() const
{
  cRow Pi(_U.size()); // just for size
  for ( size_t a = 0 ; a < _U.size() ; ++a ) 
    Pi[a] = _U[a][0]*_U[a][0];
  
  return cProbs(Pi);
}


cRow cProbModel::GetPiSqrt() const
{
  cRow PiSqrt(_U.size()); // just for size
  for ( size_t a = 0 ; a < _U.size() ; ++a ) 
    PiSqrt[a] = _U[a][0] < 0 ? -_U[a][0] : _U[a][0] ;

  return PiSqrt;
}

cRow cProbModel::GetPiSqrtRev() const
{
  cRow Pi(GetPiSqrt());
  for ( size_t i = 0 ; i < Pi.size() ; ++i )
    Pi[i] = 1.0/Pi[i];
  return (Pi);
}

void cProbModel::print(ostream& sout) const
{
  cout<<endl<<"D   ="<<endl;
  _d.Get().print(sout);
  cout<<endl<<"Phi ="<<endl;
  _Phi.print(sout);
  cout<<endl<<"U   ="<<endl;
  _U.print(sout);
  cout<<endl<<"R   ="<<endl;
  GetR().print(sout);
}

cRow cProbModel::GetPi(cSquareMatrix const& U) const
{
  assert( size() == U.size() );

  cRow Pi(U.size()); // just for size
  for ( size_t a = 0 ; a < U.size() ; ++a ) 
    Pi[a] = U[a][0]*U[a][0];

  return Pi;
}

double cProbModel::Muabt(const size_t a, const size_t b, const double t, cSquareMatrix const& U) const
{
  cRow Expt = _d.Exp(t);
  double sum = 0.;
  for ( size_t k = 0 ; k < U.size() ; ++k ) {
    sum += U[k][a]*Expt[k]*U[k][b];
  }
  return sum;
}
