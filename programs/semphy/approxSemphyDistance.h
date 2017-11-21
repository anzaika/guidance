// $Id: approxSemphyDistance.h 702 2006-05-27 17:16:58Z ninio $

#ifndef ___APPROX_SEMPHY_DISTANCE
#define ___APPROX_SEMPHY_DISTANCE

#include "computePijComponent.h"
#include "countTableComponent.h"
#include "fromCountTableComponentToDistance.h"
#include "suffStatComponent.h"
#include "semphyDistance.h"
#include "tree.h"
#include "sequenceContainer.h"


class approxSemphyDistance : public semphyDistance{

public:
	explicit approxSemphyDistance(
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const computePijGam& pij0,
				      const suffStatGlobalGam& cup,
				      const suffStatGlobalGam& cdown,
				      const VdoubleRep& cprobAtEachPos,
					  const VVdoubleRep & posteriorRateProbAtEachPos,
				      const suffStatGlobalGam& computeMarginal,
				      const Vdouble* weight,
					  const MDOUBLE toll);
	virtual ~approxSemphyDistance() {}  
	void computeDistances();
	MDOUBLE getDistance(const int nodeId1, const int nodeId2) const {
		if (nodeId1<nodeId2) return _distances[nodeId1][nodeId2];
		else return  _distances[nodeId2][nodeId1];
	}

	MDOUBLE getLikeDistance(const int nodeId1, const int nodeId2) const {
		if (nodeId1<nodeId2) return _likeDistances[nodeId1][nodeId2];
		else return  _likeDistances[nodeId2][nodeId1];
	}

	VVdouble* getDistanceTablePtr() {return &_distances;}
	VVdouble* getLikeDistanceTablePtr() {return &_likeDistances;}
	const VVdouble* getDistanceTablePtr() const {return &_distances;}
	const VVdouble* getLikeDistanceTablePtr() const {return &_likeDistances;}

	void setDistance(const int nodeId1, const int nodeId2, const MDOUBLE newVal) {
		if (nodeId1<nodeId2) _distances[nodeId1][nodeId2]=newVal;
		else _distances[nodeId2][nodeId1]=newVal;
	}
	void setLikeDistance(const int nodeId1, const int nodeId2, const MDOUBLE newVal) {
		if (nodeId1<nodeId2) _likeDistances[nodeId1][nodeId2]=newVal;
		else _likeDistances[nodeId2][nodeId1]=newVal;
	}
private:
	const suffStatGlobalGam& _cup;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	const suffStatGlobalGam& _cdown;
	const VdoubleRep& _cprobAtEachPos;
	const VVdoubleRep & _posteriorRateProbAtEachPos;

	const suffStatGlobalGam& _computeMarginal;
 	const tree& _et;
	const Vdouble* _weights;
	const MDOUBLE _toll;
	const computePijGam& _pij0;
	VVdouble _distances; // only up diagonal is used...
	VVdouble _likeDistances; // only up diagonal is used...
//	countTableComponent _ctc;
	void computeApproxDistancesNodes(
						const tree::nodeP node1,
						const tree::nodeP node2);
	void computeApproxDistancesNodesSonFather(
						const tree::nodeP nodeSon);
	void computeApproxDistancesSeparateNodes(
						const tree::nodeP nodeSon,
						const tree::nodeP nodeFather);

};

#endif

