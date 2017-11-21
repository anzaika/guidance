// cProbModel.h: interface for the cProbModel class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __CPROBMODEL_H
#define __CPROBMODEL_H

#ifdef _U
#undef _U
#endif
#ifdef _S
#undef _S
#endif

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "cEigenVals.h"
#include "cOrthonormalMatrix.h"
#include "cProbs.h"
#include "cOccuranceData.h"
#include "cSquareMatrix.h"	// Added by ClassView
#include "cAngles.h"

#include <cmath>

class cProbModel
{
public:
  cProbModel(  const cProbs& BackgroundProbs);
  virtual ~cProbModel(void);
  
  cProbs GetPi(void) const;
  cSquareMatrix GetQ(void) const;
  cSquareMatrix GetR() const;
  cRow GetPiSqrt() const;
  cRow GetPiSqrtRev() const;

  // use external U
  cRow GetPi(cSquareMatrix const& U) const;
  double Muabt(const size_t i, const size_t j, const double t, cSquareMatrix const& U) const;


  const cAngles& Angles(void) const;
  void setAngles(const cAngles& newAngles);
  void setAngle(const int k, const double newAngle);
  inline const double GetD(const int x) const {return _d.Get(x);};
  inline void SetD(const int x, const double v )  { _d.Set(x,v);};
  inline const double GetExpD(const int x, const double& t ) const {return exp(t*_d.Get(x));};
  const cEigenVals& GetD(void) const;
  const cOrthonormalMatrix& GetU(void) const;
  inline double GetU(const int a, const int b) const {return _U[a][b];};
  // 	void AdjustAllEigen(double Delta);
  // 	void RotateAll(double Delta);
  double Muabt(const size_t i, const size_t j, const double t) const;
  cSquareMatrix Mut(const double t) const;
  cSquareMatrix MutD(const double t, const cRow& D) const;
  inline size_t size() const { return _U.size(); }

  void print(ostream& sout = cout) const;
  
private:
  cEigenVals _d;
  cOrthonormalMatrix _U;
  cAngles _Phi;
  void _updateU(void);
};

#endif // __CPROBMODEL_H
