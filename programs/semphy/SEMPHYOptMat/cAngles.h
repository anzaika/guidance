#ifndef __CANGLES
#define __CANGLES
#include "cRotationMatrix.h"
#include <vector>
// cAngles is of limited responsibility
// does only what can be computed with angles only

class cAngles
{
public:
  cAngles(size_t size);
  cAngles(cProbs const& BackgroundProbs);
  cAngles(cOrthonormalMatrix const& U);

  typedef pair<size_t,size_t> tPlane;
  
  tPlane GetDim(size_t k) const;
  size_t size(void) const;
  // return the rotation matrix for angle number k.
  cRotationMatrix ROTk(const size_t k) const;
  
  // return mat*/mult_k ROT^k(\phi_k) 
  void RotateMatrix(cOrthonormalMatrix& mat) const;
  // kfirst is the first to be multiplied. if kfirst > klast: goes backwards
  void RotateMatrix(cOrthonormalMatrix& mat, int kfirst, int klast) const;

  // One of these has to go...
  void Set(const int k, const double& newangle)
  {
    _phi[k]=newangle;
  }

  void set(const int k, const double& newangle)
  {
    _phi[k]=newangle;
  }
  
  const double& get(const int k) const
  {
    return _phi[k];
  }
  
  size_t NDims(void) const;
  
  void print(ostream& sout = cout) const;
  
protected:
  void _init_phi_dim(void);
  size_t _lastK(void) const;
  size_t AngleIndex(int i, int j ) const;

private:
  size_t _size;
  vector<double> _phi;
  vector<tPlane> _dim;
};


#endif // __CANGLES

