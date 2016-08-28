#ifndef __WeightUpdateStrategy_j
#define __WeightUpdateStrategy_j

#include <vector>
using namespace std;

class tWeightUpdate;

typedef vector<double> const& Weights_cl;

class WeightUpdateStrategy_cl {
 public:
  enum WeightUpdateType { None, Random, Adversarial, Prior };

  WeightUpdateStrategy_cl( WeightUpdateType T,
			   double StartT,
			   double EndT,
			   double CoolFactor,
			   bool   UsePartialSearch )
  {
    iStart = StartT;
    iEnd = EndT;
    iCool = CoolFactor;
    iType = T;
    iPartial = UsePartialSearch;
  }
  
  tWeightUpdate* NewWeights( Weights_cl OrigWeights ) const;
  
  double NewTemp( double OldTemp ) const
  {
    double T = OldTemp*iCool;
    if( T <= iEnd )
      T = 0.0;
    return T;
  }

  bool UsePartialSearch() const
  {
    return iPartial;
  }
  
 private:
  double iStart;
  double iEnd;
  double iCool;
  bool   iPartial;
  WeightUpdateType iType;
};

#endif
