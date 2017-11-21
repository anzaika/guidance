// $Id: approxSemphyDistance.cpp 702 2006-05-27 17:16:58Z ninio $

/******************************************************************************
*                                                                             *
* This file is part of SEMPHY - Structural EM for PHYlogenetic reconstruction *
* This is unpublished proprietary source code it's authors                    *
*                                                                             *
* The contents of these coded instructions, statements and computer           *
* programs may not be disclosed to third parties, copied or duplicated in     *
* any form, in whole or in part, without the prior written permission of      *
* the authors.                                                                *
*                                                                             *
* (C) COPYRIGHT Nir Friedman, Matan Ninio,  Itsik Pe'er and Tal Pupko, 2001   *
* All Rights Reserved                                                         *
*                                                                             *
* for information please contact semphy@cs.huji.ac.il                         *
*                                                                             *
******************************************************************************/
#include "definitions.h"
#include "approxSemphyDistance.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "logFile.h"
#include <iostream>
using namespace std;

approxSemphyDistance::approxSemphyDistance(const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const computePijGam& pij0,
				      const suffStatGlobalGam& cup,
				      const suffStatGlobalGam& cdown,
				      const VdoubleRep& cprobAtEachPos,
					  const VVdoubleRep & posteriorRateProbAtEachPos,
				      const suffStatGlobalGam& computeMarginal,
				      const Vdouble* weight,
					  const MDOUBLE toll)
  : _cup(cup), _sc(sc), _sp(sp),  _cdown(cdown), _cprobAtEachPos(cprobAtEachPos),
  _posteriorRateProbAtEachPos(posteriorRateProbAtEachPos), 
		 _computeMarginal(computeMarginal),_et(et), _weights(weight) ,_toll(toll),_pij0(pij0) {
	

//	_ctc.countTableComponentAllocatePlace(_pi.alphabetSize(),_pi.categories());
	_distances.resize(_et.getNodesNum());
	_likeDistances.resize(_et.getNodesNum());
	for (int z=0; z< _et.getNodesNum(); ++z) {
		_distances[z].resize( _et.getNodesNum());
		_likeDistances[z].resize( _et.getNodesNum());
	}
}



void approxSemphyDistance::computeDistances(){
	LOG(5,<<"ENTERTING COMPUTE DISTANCES() "<<endl); 
	treeIterTopDownConst tIt(_et);
	tree::nodeP mynode = tIt.first();
	mynode = tIt.next(); // skipping the root

	for (; mynode != tIt.end(); mynode = tIt.next()) {
		computeApproxDistancesNodesSonFather(mynode);
	}

	vector<tree::nodeP> allNodes;
	_et.getAllNodes(allNodes,_et.getRoot());

	int i,j;
	for (i=0; i < allNodes.size()-1;++i) {
		for (j=i+1; j < allNodes.size(); ++j) {
			computeApproxDistancesNodes(allNodes[i],allNodes[j]);
		}
	}
}

void approxSemphyDistance::computeApproxDistancesNodes(
						const tree::nodeP node1,
						const tree::nodeP node2) {

	if (node1->father() == node2) {
		return;
	} else if (node2->father() == node1) {
		return;
	}
	else computeApproxDistancesSeparateNodes(node1,node2);
}

void approxSemphyDistance::computeApproxDistancesNodesSonFather(
						const tree::nodeP nodeSon) {

	computeCounts cc;
	countTableComponentGam ctcGam;
	cc.fillCountTableComponentGam(ctcGam,
								_sp,//const stochasticProcess& sp,
								_sc,//const sequenceContainer& sc,
								_pij0,//computePijGam& pij0,
								_cup,//const suffStatGlobalGam& cup,
								_cdown,//const suffStatGlobalGam& cdown,
								_weights,//const Vdouble * weights,
								nodeSon,//tree::nodeP nodeSon,
								_cprobAtEachPos);//const Vdouble& posProbVec
								
	const tree::nodeP nodeFather = nodeSon->father();
	const MDOUBLE tollForPairwiseDist = 0.001;
	fromCountTableComponentToDistance from1(ctcGam,_sp,tollForPairwiseDist,nodeSon->dis2father());
	from1.computeDistance();
	nodeSon->setDisToFather(from1.getDistance()); // this is actually EM step for BBL.

	if (nodeSon->id()<nodeFather->id()) {
	  _distances[nodeSon->id()][nodeFather->id()] =
	    from1.getDistance();
	  _likeDistances[nodeSon->id()][nodeFather->id()] =
	    from1.getLikeDistance();
	}
	else {
	  _distances[nodeFather->id()][nodeSon->id()] = 
	    from1.getDistance();
	  _likeDistances[nodeFather->id()][nodeSon->id()] = 
	    from1.getLikeDistance();
	}
	//	LOG(5,<<getLikeDistance(nodeSon->id(),nodeFather->id()));
	//	LOG(5,<<endl);
	
}

void approxSemphyDistance::computeApproxDistancesSeparateNodes(
						const tree::nodeP node1,
						const tree::nodeP node2) {

	countTableComponentGam ctcGam;
	ctcGam.countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories());

	if (_weights!=NULL) {
		for (int pos=0; pos< _sc.seqLen(); ++pos) {
			if ((*_weights)[pos] == 0.0) continue; 
			for (int alph1 =0; alph1< _sp.alphabetSize(); ++alph1) {
				for (int alph2 =0; alph2< _sp.alphabetSize(); ++alph2) {
					for (int rateCategor =0; rateCategor< _sp.categories(); ++rateCategor) {
						doubleRep tmp = 
						  _computeMarginal.get(
										 pos,
										 rateCategor,
										 node1->id(),
										 alph1
										 ) *
						  _computeMarginal.get(
										 pos,
										 rateCategor,
										 node2->id(),
										 alph2
										 )
										  * _posteriorRateProbAtEachPos[pos][rateCategor];
											
						ctcGam.addToCounts(alph1,alph2,rateCategor,convert(tmp)*(*_weights)[pos]);
					}
				}
			}
		}
	}else{		// _weights==NULL
		for (int pos=0; pos< _sc.seqLen(); ++pos) {
			for (int alph1 =0; alph1< _sp.alphabetSize(); ++alph1) {
				for (int alph2 =0; alph2< _sp.alphabetSize(); ++alph2) {
					for (int rateCategor =0; rateCategor< _sp.categories(); ++rateCategor) {
						doubleRep tmp = 
						  _computeMarginal.get(pos,rateCategor,node1->id(),alph1)
										 *
						  _computeMarginal.get(pos,rateCategor,node2->id(),alph2)
						  * 
						 _posteriorRateProbAtEachPos[pos][rateCategor];

						  ctcGam.addToCounts(alph1,alph2,rateCategor,convert(tmp));
					}
				}
			}
		}
	}
	MDOUBLE initialGuess = getDistance(node1->id(),node2->id());
	if (initialGuess==0) initialGuess=0.031;
	const MDOUBLE tollForPairwiseDist = 0.001;
	fromCountTableComponentToDistance from1(ctcGam,_sp,tollForPairwiseDist,initialGuess);
	from1.computeDistance();
	if (node1->id()<node2->id()) {
		_distances[node1->id()][node2->id()] = 
			from1.getDistance();
		_likeDistances[node1->id()][node2->id()] = 
			from1.getLikeDistance();
	}
	else {
		_distances[node2->id()][node1->id()] = 
			from1.getDistance();
		_likeDistances[node2->id()][node1->id()] = 
			from1.getLikeDistance();
	}
}


