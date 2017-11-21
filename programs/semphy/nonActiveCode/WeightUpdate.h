#ifndef __WeightUpdate_h
#define __WeightUpdate_h

//#include <spin.h>
//#include <MainLib.h> //CHANGE
#include <vector> //change
using namespace std; //change
/** This class is an annealing weight update interface. It keeps record
    of the original weights and previous weights and returns a new
    weight vector given a vector of gradients for the weights. 
    Update is governed by the temperature
*/
class tWeightUpdate {

 public:

  /// default constructor
  tWeightUpdate(const vector<double>& origWeights) : 
    _origWeights(origWeights), _prevWeights(origWeights), _temp(0.0) {};

  /// default destructor
  virtual ~tWeightUpdate() {};

  /// update the weights using gradients
  virtual vector<double> const& UpdateWeights(double newtemp,vector<double> const& gradient) = 0;

  /// return the last weight vector
  vector<double> const& CurrentWeights() { return _prevWeights; };

  /// return original weight vector
  vector<double> const& OrigWeights() { return _origWeights; };
  ///
  double GetTemp() { return _temp; };

  /// set the temperature
  void SetTemperature(double temp);
  
 protected:
  /// get non const reference to prev weights
  vector<double>& PrevWeights() { return _prevWeights; };
  
 private:
  vector<double> _origWeights;
  vector<double> _prevWeights;
  double _temp;
};

/** This class implement the interface of tWeightUpdate by
    random updating of weights centered around the original weights
    in proportion to the temperature.
 */
class tWeightUpdateRandom : public tWeightUpdate {
 public:
  /// default constructor
  tWeightUpdateRandom(const vector<double>& origWeights) : tWeightUpdate(origWeights) {};

  /// default destructor
  virtual ~tWeightUpdateRandom() {};

  /// update the weights ignoring gradients
  virtual vector<double> const& UpdateWeights(double newtemp,vector<double> const& gradient);
 private:
};

/** This class implement the interface of tWeightUpdate by
    doing a gradient update that is governed by the temperature
    and damping factors for the original and previous weights. 
    The learning rate determines the magnitude of update.
 */
class tWeightUpdateGradient : public tWeightUpdate {
 public:
  /// default constructor
  tWeightUpdateGradient(const vector<double>& origWeights) : 
    tWeightUpdate(origWeights), _lRate(1.0), _origDamping(1.0), _prevDamping(1.0) {};
  /// default destructor
  virtual ~tWeightUpdateGradient() {};
  /// update the weights using gradients
  virtual vector<double> const& UpdateWeights(double newtemp,vector<double> const& gradient);
  
  /// set the learning rate for multiplicative update
  void SetLearningRate(double lRate) { _lRate = lRate; };
  /// set damping factor with respect to original weights
  void SetOrigDamping(double origDamping) { _origDamping = origDamping; };
  /// set damping factor with respect to previous weights
  void SetPrevDamping(double prevDamping) { _prevDamping = prevDamping; };
 private:
  double _lRate;
  double _origDamping;
  double _prevDamping;
};

#endif




