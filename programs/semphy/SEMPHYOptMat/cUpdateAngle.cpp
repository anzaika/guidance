// cUpdateAngle.cpp: implementation of the cUpdateAngle class.
//
//////////////////////////////////////////////////////////////////////
#include <map>
#include "cCoeffMatrices.h"
#include "cUpdateAngle.h"
#include <cmath>

using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cUpdateAngle::cUpdateAngle(cDataProbModel& DataProbModel):
	_data_prob_model(DataProbModel), _deriv_uxy(_data_prob_model. DerivAllUxy())
{
}

cUpdateAngle::~cUpdateAngle()
{

}

double
cUpdateAngle::ComputedPhi(size_t k)
{
  cCoeffMatrices CoeffMatrices( _data_prob_model.ProbModel().Angles(),k );
  map<cDRotationMatrix::tDrotType,double> Theta;

  //  cout<<"Theta=[";

  for ( cDRotationMatrix::tDrotType RotType = (cDRotationMatrix::tDrotType)0 
	  ; RotType < cDRotationMatrix::NROTT 
	  ; RotType=(cDRotationMatrix::tDrotType) ++((int) RotType)) 
  {
    cSquareMatrix MMat(CoeffMatrices.CoeffMatrix(RotType));
    Theta[RotType] = _deriv_uxy.DotProduct(MMat);

  }
  double phi = _data_prob_model.ProbModel().Angles().get(k);
  return Theta[cDRotationMatrix::SIN]*cos( phi ) - Theta[cDRotationMatrix::COS]*sin( phi );
}

double cUpdateAngle::ComputeNewPhi(size_t k)
{	 
  cCoeffMatrices CoeffMatrices( _data_prob_model.ProbModel().Angles(),k );
  map<cDRotationMatrix::tDrotType,double> Theta;

  //  cout<<"Theta=[";

  for ( cDRotationMatrix::tDrotType RotType = (cDRotationMatrix::tDrotType)0 
	  ; RotType < cDRotationMatrix::NROTT 
	  ; RotType=(cDRotationMatrix::tDrotType) ++((int) RotType)) {
    cSquareMatrix MMat(CoeffMatrices.CoeffMatrix(RotType));

    Theta[RotType] = _deriv_uxy.DotProduct(MMat);
    //    cout<<Theta[RotType]<<" ";
  }
  //  cout<<"]"<<endl;

  // Sanity check
  //  cout<<"deriv_uxy"<<endl;
  //  _deriv_uxy.print(cout);

  {
    double phi = _data_prob_model.ProbModel().Angles().get(k);
    double dPhi = Theta[cDRotationMatrix::SIN]*cos( phi ) - Theta[cDRotationMatrix::COS]*sin( phi );
    cout << "dPhi = " << dPhi << "\n";
    // Print empirical dPhi

  }
  return atan(Theta[cDRotationMatrix::SIN]/Theta[cDRotationMatrix::COS]);
}

cAngles cUpdateAngle::ComputeNewPhi()
{
	cAngles Result(_data_prob_model.ProbModel().Angles().size());
	for ( size_t k = 0 ; k < Result.NDims() ; ++k ) {
		Result.Set(k,ComputeNewPhi(k));
	}
	return Result;
}


cAngles cUpdateAngle::ComputeNewPhiNoPiChange()
{
	cAngles Result(_data_prob_model.ProbModel().Angles().size());
	for ( size_t k = Result.size() -1; k < Result.NDims() ; ++k ) {
		Result.Set(k,ComputeNewPhi(k));
	}
	return Result;
}
