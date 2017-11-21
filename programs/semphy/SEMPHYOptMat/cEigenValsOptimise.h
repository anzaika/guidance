#ifndef   __CEIGENVALSOPTIMISE
#define   __CEIGENVALSOPTIMISE


#include "cDataProbModel.h"

#define MAXEIGENVAL (-0.01)
#define MINEIGENVAL (-0.6)

class cEigenValsOptimise {
public:
  //  double optimize(const int x) const;
  double optimize(const int x, const double tol) const;

  cEigenValsOptimise(const cDataProbModel& DM): _DM(DM){};
  
  double LL_dx(const double& dx, size_t x) const;
  double LL_dx(const cRow& D) const;
  double dLL_ddx(const cRow& D, const int x) const;

private:


  const cDataProbModel& _DM;

};




#endif // __CEIGENVALSOPTIMISE
