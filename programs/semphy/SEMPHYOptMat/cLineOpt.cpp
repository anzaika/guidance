// cLineOpt.cpp: implementation of the cLineOpt class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "cLineOpt.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cLineOpt::cLineOpt(const cOccuranceData& OccuranceData,
				   const cProbs& BackgroundProbs) :
	_prob_model(OccuranceData, BackgroundProbs)
{

}

cLineOpt::~cLineOpt()
{

}

#define MIN_LL_IMPROVEMENT .0001
#define INIT_DELTA_EIG .01
#define INIT_DELTA_ROT .01
#define EIG_OPT_INC_FACTOR 1.2
#define ROT_OPT_INC_FACTOR 1.2
#define EIG_OPT_DEC_FACTOR 5
#define ROT_OPT_DEC_FACTOR 5

void cLineOpt::Opitmize(void)
{
	double DeltaEig = INIT_DELTA_EIG;
	double DeltaRot = INIT_DELTA_ROT;
	double NewLL = _prob_model.NonConstLL();
	double PrevLL;
	do {
		PrevLL = NewLL;
		DeltaEig = ImproveEig(DeltaEig);
		DeltaRot = ImproveRot(DeltaRot);
		NewLL = _prob_model.NonConstLL() - PrevLL;		
	} while ( NewLL - PrevLL > MIN_LL_IMPROVEMENT );
}

double cLineOpt::ImproveEig(double Delta)
{
	double NewLL = _prob_model.NonConstLL();
	double PrevLL;
	do {
		PrevLL = NewLL;
		_prob_model.AdjustAllEigen(Delta);
		NewLL = _prob_model.NonConstLL();		
		Delta *= EIG_OPT_INC_FACTOR;
	} while ( NewLL > PrevLL );
	// last move spoilt likelihood. Correct it!
	Delta /= EIG_OPT_INC_FACTOR;
	_prob_model.AdjustAllEigen(-Delta);
	return Delta/EIG_OPT_DEC_FACTOR;
}

double cLineOpt::ImproveRot(double Delta)
{
	double NewLL = _prob_model.NonConstLL();
	double PrevLL;
	do {
		PrevLL = NewLL;
		_prob_model.RotateAll(Delta);
		NewLL = _prob_model.NonConstLL();		
		Delta *= ROT_OPT_INC_FACTOR;
	} while ( NewLL > PrevLL );
	// last move spoilt likelihood. Correct it!
	Delta /= ROT_OPT_INC_FACTOR;
	_prob_model.RotateAll(-Delta);
	return Delta/ROT_OPT_DEC_FACTOR;
}

	
const cEigenVals& cLineOpt::Getd(void) const
{
	return _prob_model.Getd();
}
	
const cOrthonormalMatrix& cLineOpt::GetU(void) const
{
	return _prob_model.GetU();
}

