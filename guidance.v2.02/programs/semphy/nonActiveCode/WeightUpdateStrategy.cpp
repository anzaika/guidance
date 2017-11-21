#include "WeightUpdate.h" //change
#include "WeightUpdateStrategy.h" //change
#include <cassert> //change
tWeightUpdate*
WeightUpdateStrategy_cl::NewWeights( Weights_cl OrigWeights ) const
{
  tWeightUpdate* W;
  switch( iType )
  {
  case None:
  case Random:
  case Prior:
    W = new tWeightUpdateRandom( OrigWeights );
    break;
  case Adversarial:
    W = new tWeightUpdateGradient( OrigWeights );
    break;
  default:
    assert(false);
  }
  W->SetTemperature( iStart/iCool );
  if( iType == None || iType == Prior )
    W->SetTemperature( 0.0 );
    
  return W;
}
