// 	$Id: constantSemphyDistance.h 922 2006-09-21 09:36:33Z ninio $	

#ifndef ___CONSTANT_SEMPHY_DISTANCE
#define ___CONSTANT_SEMPHY_DISTANCE

#include "semphyDistance.h"

class constantSemphyDistance : public semphyDistance{
public:
  explicit constantSemphyDistance( const MDOUBLE & dist = 1.0, const MDOUBLE & like = 0.0):_dist(dist),_like(like) {};
  virtual ~constantSemphyDistance() {} ; 
  void computeDistances(){};
  MDOUBLE getDistance(const int nodeId1, const int nodeId2) const {return _dist;}
  MDOUBLE getLikeDistance(const int nodeId1, const int nodeId2) const { return _like; }
  VVdouble* getDistanceTablePtr() {return NULL;}
  VVdouble* getLikeDistanceTablePtr() {return NULL;}
  const VVdouble* getDistanceTablePtr() const {return NULL;}
  const VVdouble* getLikeDistanceTablePtr() const {return NULL;}

  void setDistance(const int nodeId1, const int nodeId2, const MDOUBLE newVal) { _dist=newVal;}
  void setLikeDistance(const int nodeId1, const int nodeId2, const MDOUBLE newVal) { _like=newVal;}
private:
  MDOUBLE _dist, _like;
};

#endif

