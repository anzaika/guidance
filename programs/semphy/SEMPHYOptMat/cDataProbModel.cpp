// cDataProbModel.cpp: implementation of the cDataProbModel class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cDataProbModel.h"
#include "cUpdateAngle.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cDataProbModel::cDataProbModel(	const cOccuranceData& OccuranceData,
							   cProbModel& ProbModel)
: _prob_model(ProbModel), _D(OccuranceData)
{

}

cDataProbModel::~cDataProbModel()
{

}



double cDataProbModel::NonConstLL(cSquareMatrix const& U) const
{
  double sum = 0.;
  for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) 
  {
    double t = _D.Time(Which_t);
    //    cout <<"t="<<t<<endl;
    cSquareMatrix const& CurrentMat(_D.Matrix(Which_t));

    for ( size_t i = 0 ; i < size() ; ++i ) {
      for ( size_t j = 0 ; j < size() ; ++j ) {
	sum += CurrentMat[i][j]*log(ProbModel().Muabt(i,j,t,U) );
      }
    }
  }
  return sum;
}

double cDataProbModel::LL(cSquareMatrix const& U) const
{
  double sum = NonConstLL(U);
  cRow pi(ProbModel().GetPi(U));
  cRow logPi(pi.size());
  for (int i=0;i<pi.size();++i)
    logPi[i]=log(pi[i]);		// log pi
  
  // inefficient implementation!
  for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) {
    const cSquareMatrix CurrentMat(_D.Matrix(Which_t));
    for ( size_t a = 0 ; a < size() ; ++a ) {
      for ( size_t b = 0 ; b < size() ; ++b ) {
	sum += CurrentMat[a][b]*0.5*(logPi[b] - logPi[a]);
      }
    }
  }
  return sum;
}


double cDataProbModel::NonConstLL() const
{
	double sum = 0.;
	for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) {
		const cSquareMatrix CurrentMat(_D.Matrix(Which_t));
		for ( size_t i = 0 ; i < size() ; ++i ) {
			for ( size_t j = 0 ; j < size() ; ++j ) {
				sum += CurrentMat[i][j]*log(MuabWhicht(i,j,Which_t));
			}
		}
	}
	return sum;
}

double cDataProbModel::LL(void) const
{
  double sum = NonConstLL();
  cRow pi(ProbModel().GetPi());
  cRow logPi(pi.size());
  for (int i=0;i<pi.size();++i)
    logPi[i]=log(pi[i]);		// log pi
  
  // inefficient implementation!
  for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) {
    const cSquareMatrix CurrentMat(_D.Matrix(Which_t));
    for ( size_t a = 0 ; a < size() ; ++a ) {
      for ( size_t b = 0 ; b < size() ; ++b ) {
	sum += CurrentMat[a][b]*0.5*(logPi[b] - logPi[a]);
      }
    }
  }
  return sum;
}

double cDataProbModel::DerivUxy(size_t x, size_t y) const
{
  cRow pi(ProbModel().GetPi());
	double sum = 0.;
	double dx = _prob_model.GetD(x);
	for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) {
	  double ExpVal = _D.Time(Which_t)*dx;
	  double etdx = exp(ExpVal);
	  const cSquareMatrix& CurrentMat(_D.Matrix(Which_t));
	  for ( size_t a = 0 ; a < size() ; ++a ) {
	    //	    	    cout << ProbModel().GetU(x,a)<<"*"<<(CurrentMat[a][y]+CurrentMat[y][a])<<"*"<<etdx<<"/"<<MuabWhicht(x,y,Which_t);
		    double tmp = ProbModel().GetU(x,a)*(CurrentMat[a][y]+CurrentMat[y][a])*etdx
	      /MuabWhicht(y,a,Which_t);
		    //		    cout <<"="<<tmp<<endl;
		      sum += tmp;
	  }
	}
	if (y == 0)
	{
	  sum += _D.getSstar(x)/(ProbModel().GetU(x,0));
	}
	return sum;
}

double cDataProbModel::DerivRotxy(size_t i, size_t j) const
{
	double sum = 0.;
	for ( size_t k = 0 ; k < size() ; ++k ) {
		sum += DerivUxy(i,k);
		sum -= DerivUxy(j,k);
	}
	return sum;
}

cRow
cDataProbModel::DerivAllEigen() const
{
  cRow Res( size() );
  for( size_t a = 0; a < size(); a++ )
    Res[a] = DerivEigen(a);
  return Res;
}
  

double cDataProbModel::DerivEigen(size_t x) const
{
	double sum = 0.;
	for ( size_t Which_t = 0 ; Which_t < _D.size() ; ++Which_t ) {
		const cSquareMatrix CurrentMat(_D.Matrix(Which_t));
		double t = _D.Time(Which_t);
		double etdx = exp(t*ProbModel().GetD(x));
		for ( size_t a = 0 ; a < size() ; ++a ) {
		  for ( size_t b = 0 ; b < size() ; ++b ) {
		    sum += t*etdx*CurrentMat[a][b]*ProbModel().GetU(x,a)*ProbModel().GetU(x,b)
		      /_prob_model.Muabt(a,b,t);
		  }
		}
	}
	return sum;
}


cSquareMatrix cDataProbModel::DerivAllUxy(void) const
{
	cSquareMatrix Result(size());
 	for ( size_t i = 0 ; i < size() ; ++i ) {
	 	for ( size_t j = 0 ; j < size() ; ++j ) {
			Result.Set(i,j,DerivUxy(i,j));
 		}
	}
	return Result;
}

#ifdef JUNK
cSquareMatrix cDataProbModel::DerivAllRotxy(void) const
{
	cSquareMatrix Result(size());
 	for ( size_t i = 0 ; i < size() ; ++i ) {
	 	for ( size_t j = 0 ; j < size() ; ++j ) {
			Result.Set(i,j,DerivRotxy(i,j));
 		}
	}
	return Result;
}
#endif


cRow 
cDataProbModel::DerivAllAngles(void)
{
  cUpdateAngle UA( *this );

  cRow R( ProbModel().Angles().NDims() );
  for( size_t k = 0; k < ProbModel().Angles().NDims(); k++ )
    R[k] = UA.ComputedPhi(k);

  return R;
}

cProbModel& cDataProbModel::ProbModel()
{
	return _prob_model;
}

cProbModel const& cDataProbModel::ProbModel() const
{
	return _prob_model;
}

size_t cDataProbModel::size() const
{
	return ProbModel().GetU().size();
}

double cDataProbModel::MuabWhicht(const size_t a, 
								  const size_t b, 
								  const size_t Which_t) const
{
	return ProbModel().Muabt(a,b,_D.Time(Which_t));
}
