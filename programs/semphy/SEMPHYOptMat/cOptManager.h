// cOptManager.h: interface for the cOptManager class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COPTMANAGER_H__9FF6EA37_611D_485E_81B0_18D12B2D4816__INCLUDED_)
#define AFX_COPTMANAGER_H__9FF6EA37_611D_485E_81B0_18D12B2D4816__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cProbs.h"
#include "cDataProbModel.h"


class cOptManager{
public:

  typedef enum {PARALLEL, SEQUENTIAL} tModeOpt;

  struct tOptParams {
    bool      FreezeBackground;
    tModeOpt  OptType;
    double    StopCriteria;
    int       MaxIter;

    tOptParams() {
      FreezeBackground = true;
      OptType = SEQUENTIAL;
      StopCriteria = 0.00001;
      MaxIter = 100;
    }
  };
  
  cOptManager( cOccuranceData const& OccuranceData, 
	       cProbs& BackgroundProbs,
	       tOptParams Params = tOptParams() );

  cOptManager( cOccuranceData const& OccuranceData, 
	       cProbModel const& Model,
	       tOptParams Params = tOptParams() );

  virtual ~cOptManager();
  
  virtual
  cProbModel operator()(void);
  
protected:
  
  size_t First_Angle(void) const;

  void IterationAngleSequential(cDataProbModel& DM);

  void IterationAngleParallel(cDataProbModel& DM);

  void IterationEigen(cDataProbModel& DM);

  
private:
  
  // disable copy & assignment
  cOptManager(const cOptManager&);
  const cOptManager& operator=(const cOptManager&);

  //
  tOptParams _params;
  cOccuranceData const& _data;
  cProbModel _initModel;
};

#endif // !defined(AFX_COPTMANAGER_H__9FF6EA37_611D_485E_81B0_18D12B2D4816__INCLUDED_)
