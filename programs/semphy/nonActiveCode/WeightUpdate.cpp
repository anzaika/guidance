//#include <MathFunctions.h>
#include <cassert> //CHANGE
#include "WeightUpdate.h" //CHANGE
//#include <RandomProb.h>
#include "RandomGenerator.h" //CHANGE

extern tRandomGenerator _RandomProbGenerator;

const double TINY = 10e-7;

void tWeightUpdate::SetTemperature(double temp) 
{ 
  _temp = temp; 
  if ( _temp == 0.0 ) { // need to restore to original weights
    for ( int i=0 ; i<_origWeights.size() ; i++ )
      _prevWeights[i] = _origWeights[i];
  }
};

vector<double> const& tWeightUpdateRandom::UpdateWeights(double newtemp,vector<double> const& gradient)
{
  //  cerr << "Updating weight for random\n";
  SetTemperature(newtemp);
  if ( newtemp == 0.0 )
    return CurrentWeights();

  vector<double>& prevWeights = PrevWeights();
  vector<double> const& origWeights = OrigWeights();
  //mn130103  double scale = sqrt(newtemp);
  double var = 1.0/newtemp;
  double wsum = 0.0;
  double origsum = 0.0;
  int i;
  // go over the instances and update new weights
  for ( i=0 ; i<prevWeights.size() ; i++ ) {
    //mn130103    double noise = 0.0;
    double weight = 0.0;
    // sampling is according to original weight
    // This assumes original weights are whole numbers!!!
    for ( int w=1 ; w <= origWeights[i] ; w++ )
	weight += _RandomProbGenerator.SampleGamma(var) * newtemp;
    if ( weight < TINY )
      weight = TINY;
    prevWeights[i] = weight;
    origsum += origWeights[i];
    wsum += weight;
  }
  // normalize
  double fact = origsum / wsum;
  for ( i=0 ; i<prevWeights.size() ; i++ ) {
    prevWeights[i] *= fact;
    //    cerr << "Weight " << i << "=" << prevWeights[i] << endl;
  }

  return CurrentWeights();
}

vector<double> const& tWeightUpdateGradient::UpdateWeights(double newtemp,vector<double> const& gradient)
{
  //  cerr << "Updating weight for gradient\n";
  // make sure we got all the gradients
  assert(gradient.size()==OrigWeights().size());

  cerr << "new temp = "<<newtemp<<endl;
  SetTemperature(newtemp);
  if ( newtemp == 0.0 )
    return CurrentWeights();

  vector<double>& prevWeights = PrevWeights();
  vector<double> const& origWeights = OrigWeights();
  double beta = _origDamping/newtemp;
  double gamma = _prevDamping/newtemp;
  double sfact = 1.0/(beta+gamma);
  double bfact = beta*sfact;
  double gfact = gamma*sfact;
  double origsum = 0.0;
  double maxw = -HUGE_VAL;
  // go over the instances and update new weights
  int i;
  for ( i=0 ; i<prevWeights.size() ; i++ ) {
    origsum += origWeights[i];
    double weight = bfact*log(origWeights[i]) + gfact*log(prevWeights[i]) - sfact*_lRate*gradient[i];
    prevWeights[i] = weight;
    if ( weight > maxw )
	  maxw = weight;
  }
  // normalize
  double wsum = 0.0;
  for ( i=0 ; i<prevWeights.size() ; i++ ) {
    prevWeights[i] -= maxw;
    prevWeights[i] = exp(prevWeights[i]);
    if ( prevWeights[i] < TINY )
      prevWeights[i] = TINY;
    wsum += prevWeights[i];
  }
  double f = origsum / wsum;
  for ( i=0 ; i<prevWeights.size() ; i++ ) {
    prevWeights[i] *= f;
    cerr << "Weight " << i << "=" << prevWeights[i] << "("<<gradient[i]<<")"<<endl;
  }

  return CurrentWeights();
}
