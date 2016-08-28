// $Id: exactSemphyDistance.h 409 2005-06-28 13:12:24Z ninio $

#ifndef ___EXACT_SEMPHY_DISTANCE
#define ___EXACT_SEMPHY_DISTANCE
/*
#include "stdAfxLocal.h"
#include "stdAfx.h"
#include "fromCountTableComponentToDistance.h"
#include "suffStatComponent.h"
#include "computeProbOfEachPos.h"
#include "semphyDistance.h"
#include "positionInfo.h"

class exactSemphyDistance : public semphyDistance{

public:
	explicit exactSemphyDistance::exactSemphyDistance(const tree& et,
							  const positionInfo &pi,
							  const suffStatComponent& cup,
							  const suffStatComponent& cdown,
							  const computeProbOfEachPos& cprobAtEachPos,
							  const suffStatComponent& computeMarginal,
							  const Vdouble* weight,
							  const MDOUBLE toll);
	virtual ~exactSemphyDistance();
	void computeDistances() ;
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
	const suffStatComponent& _cup;
	const positionInfo& _pi;
	const suffStatComponent& _cdown;
	const computeProbOfEachPos& _cprobAtEachPos;
	const suffStatComponent& _computeMarginal;
 	const tree& _et;
	const Vdouble* _weight;
	const Vdouble* _tmpWeight;
	const MDOUBLE _toll;

	VVdouble _distances; // only up diagonal is used...
	VVdouble _likeDistances; // only up diagonal is used...


	vector < vector < countTableComponent *> > _vvctcp;
	vector < vector < countTableComponent *> > _tmpvvctcp;
  

	void computeExactDistancesNodes(const tree::nodeP node1,
						const tree::nodeP node2);
	void computeExactDistancesNodesSonFather(const tree::nodeP nodeSon,
						 const int pos);
	void computeExactDistancesSeparateNodes(const tree::nodeP newNode,
						const tree::nodeP oldNode);

	void acculate_ctc(const int from, const int to);
	void exactSemphyDistance::ctcFromTwoCtc(countTableComponent* resCtc, const countTableComponent& newCtc, const countTableComponent& oldCtc);
	void exactSemphyDistance::ctcFromTwoCtcReverse(countTableComponent* resCtc, const countTableComponent& newCtc, const countTableComponent& oldCtc);
};
*/
#endif

