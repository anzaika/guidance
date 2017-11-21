// cUpdateAngle.h: interface for the cUpdateAngle class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CUPDATEANGLE_H__2CAC1A30_D49C_4862_9DB4_B41021FC26DB__INCLUDED_)
#define AFX_CUPDATEANGLE_H__2CAC1A30_D49C_4862_9DB4_B41021FC26DB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cEigenVals.h"
#include "cDataProbModel.h"
#include "cAngles.h"	// Added by ClassView

class cUpdateAngle  
{
public:
  cUpdateAngle(cDataProbModel& ProbModel);

  virtual ~cUpdateAngle();

  cAngles ComputeNewPhi(void);

  cAngles ComputeNewPhiNoPiChange(void);

  double ComputedPhi(size_t k);

  double ComputeNewPhi(size_t k);


private:
  cDataProbModel& _data_prob_model;

  cSquareMatrix _deriv_uxy;
};

#endif // !defined(AFX_CUPDATEANGLE_H__2CAC1A30_D49C_4862_9DB4_B41021FC26DB__INCLUDED_)
