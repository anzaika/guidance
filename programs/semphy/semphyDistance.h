// $Id: semphyDistance.h 409 2005-06-28 13:12:24Z ninio $

#ifndef ___SEMPHY_DISTANCE
#define ___SEMPHY_DISTANCE

#include "definitions.h"


class semphyDistance{

public:

	virtual void computeDistances() =0;
	virtual ~semphyDistance() =0;

	virtual VVdouble* getDistanceTablePtr()=0;
	virtual const VVdouble* getDistanceTablePtr() const =0;
	virtual VVdouble* getLikeDistanceTablePtr()=0;
	virtual MDOUBLE getDistance(const int nodeId1, const int nodeId2) const =0;
	virtual MDOUBLE getLikeDistance(const int nodeId1, const int nodeId2) const =0;
};

#endif
