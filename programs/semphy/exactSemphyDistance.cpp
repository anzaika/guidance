// $Id: exactSemphyDistance.cpp 409 2005-06-28 13:12:24Z ninio $

#include "exactSemphyDistance.h"
/*
exactSemphyDistance::exactSemphyDistance(const tree& et,
					 const positionInfo &pi,
					 const suffStatComponent& cup,
					 const suffStatComponent& cdown,
					 const computeProbOfEachPos& cprobAtEachPos,
					 const suffStatComponent& computeMarginal,
					 const Vdouble* weight,
					 const MDOUBLE toll):
	_cup(cup), _pi(pi),
	_cdown(cdown), _cprobAtEachPos(cprobAtEachPos),
	_computeMarginal(computeMarginal), _et(et), _weight(weight), _tmpWeight(NULL),_toll(toll){
	_vvctcp.resize(_et.iNodes());
	int i=0;
	for (i=0;i<_et.iNodes();++i) _vvctcp[i].resize(_et.iNodes(),NULL);
	
	_tmpvvctcp.resize(_et.iNodes());
	for (i=0;i<_et.iNodes();++i) _tmpvvctcp[i].resize(_et.iNodes(),NULL);
	
	VVresize(_distances,_et.iNodes(),_et.iNodes());
	VVresize(_likeDistances,_et.iNodes(),_et.iNodes());
	if (weight==NULL) {
	  _tmpWeight = new Vdouble(_pi.seqLen(),1.0);
	  _weight = _tmpWeight;
	}
}

exactSemphyDistance::~exactSemphyDistance(){
	for (int i=0;i<_vvctcp.size();++i){
		for (int j=0;j<_vvctcp[i].size();++j){
			if (_vvctcp[i][j]!=NULL) delete  _vvctcp[i][j];
			if (_tmpvvctcp[i][j]!=NULL)	delete  _tmpvvctcp[i][j];
		}
	}
	if (_tmpWeight != NULL) delete _tmpWeight;
}

void exactSemphyDistance::computeDistances(){
	vector<tree::nodeP> allNodes;
	_et.getAllNodesBFS(allNodes,_et.getRoot()); 
	for (int pos=0; pos< _pi.seqLen(); ++pos) {
		if (_weight!=NULL && (*_weight)[pos] == 0) continue;	  
		for (vector<tree::nodeP>::iterator newNode=allNodes.begin(); newNode != allNodes.end();++newNode) {
			if (_et.isRoot(*newNode)==false) {
				computeExactDistancesNodesSonFather(*newNode,pos);
				for (vector<tree::nodeP>::iterator oldNode=allNodes.begin(); oldNode != newNode;++oldNode) {
					if ((*newNode)->father() != *oldNode) {// don't calculate from you to your father, just becouse your father was allready seen
						computeExactDistancesSeparateNodes(*newNode, *oldNode); // no need to specify position, as we do not use it explicitly
					}
				}
			}
		}
	}
	// we now calculated the CTCs, we still need to compute the distances

	for (vector<tree::nodeP>::iterator fromNode=allNodes.begin(); fromNode != allNodes.end();++fromNode) {
		for (vector<tree::nodeP>::iterator toNode=allNodes.begin(); toNode != fromNode;++toNode) {
			fromCountTableComponentToDistance from1(*_vvctcp[(*fromNode)->id()][(*toNode)->id()],
													*_pi.stocProcess(0),
													_toll,
													(*fromNode)->dis2father());
			from1.computeDistance();
			if ((*fromNode)->id()<(*toNode)->id()){
				_distances[(*fromNode)->id()][(*toNode)->id()]=from1.getDistance();
				_likeDistances[(*fromNode)->id()][(*toNode)->id()]=from1.getLikeDistance();
			} else {
				_distances[(*toNode)->id()][(*fromNode)->id()] = from1.getDistance();
				_likeDistances[(*toNode)->id()][(*fromNode)->id()]=from1.getLikeDistance();
			}
		}
	}
}

void exactSemphyDistance::acculate_ctc(const int from, const int to){
	if (_vvctcp[from][to]==NULL){
		_vvctcp[from][to]=new countTableComponent;
		_vvctcp[from][to]->countTableComponentAllocatePlace(_pi.alphabetSize(),_pi.categories());
		_vvctcp[from][to]->zero();
	}
	for (int i=0;i<_vvctcp[from][to]->_countValues.size();++i){
		for (int j=0;j<_vvctcp[from][to]->_countValues[0].size();++j){
			for (int k=0;k<_vvctcp[from][to]->_countValues[0][0].size();++k){
				_vvctcp[from][to]->_countValues[i][j][k] += _tmpvvctcp[from][to]->_countValues[i][j][k];
			}
		}
	}
}



void exactSemphyDistance::computeExactDistancesNodesSonFather(
						const tree::nodeP nodeSon,
						int pos) {
	const tree::nodeP& nodeFather = nodeSon->father();

	int from=nodeSon->id();
	int to=nodeFather->id();
	if (_tmpvvctcp[from][to]==NULL){
		_tmpvvctcp[from][to]=new countTableComponent;
		_tmpvvctcp[from][to]->countTableComponentAllocatePlace(_pi.alphabetSize(),_pi.categories());
		_tmpvvctcp[from][to]->zero();
	}
	
	_tmpvvctcp[from][to]->zero(); // clear the right ctc
	for (int alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
		for (int alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
			for (int rate =0; rate< _pi.categories(); ++rate) {
				MDOUBLE tmp = _cup.get(from,pos,rate,alph1) *
				_cdown.get(from,pos,rate,alph2) *
				_pi.pij(0)->getPij(from,alph1,alph2,rate)*
				_pi.stocProcess(0)->freq(alph1)/
				_cprobAtEachPos.getProb(pos);
		  
				_tmpvvctcp[from][to]->addToCounts(alph1,alph2,rate,tmp*(*_weight)[pos]);
			}
		}
	}
	//_tmpvvctcp[from][to]->printTable();
	//exit(4);
	acculate_ctc(from,to);
}
// BUGHERE???
// this function will multiply two CTC's to form a new one.
  // I'm doing this while assuming that we allways keep the ctc's from the newer nodes to the older ones.
  // this is consistant with computing from node to father, but not the reverse.
  // doing it the other way around whould make more sence computational, as we need to computer the reverse of one 
  // of the ctc's, namely the secound one.  as we allways computer the pairs while using the new CTC, 
  // we can reverse it once and save a lot of work.  not implinted here.

// MATAN VERSION...
void exactSemphyDistance::ctcFromTwoCtc(countTableComponent* resCtc,
										const countTableComponent& newCtc,
										const countTableComponent& oldCtc){
	countTableComponent tmp(oldCtc);
	//  now, we sum each line to computer p(b)=SIGMA_c (P(b,c)) so we can compute P(c|b)=p(b,c)/p(b)
	int rateCategor =0;
	for (rateCategor =0; rateCategor< _pi.categories(); ++rateCategor) {
		for (int alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
			MDOUBLE sum=0.0;
			int alph2=0;
			for (alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
				sum += tmp._countValues[alph1][alph2][rateCategor];
			}
			for (alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
				tmp.setCount(alph1, alph2, rateCategor, oldCtc._countValues[alph1][alph2][rateCategor]/(sum+EPSILON));
			}
		}
	}
	for (rateCategor =0; rateCategor< _pi.categories(); ++rateCategor) {
		for (int alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
			for (int alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
				MDOUBLE sum_location=0.0; 
				for (int alph_mid = 0; alph_mid < _pi.alphabetSize() ; ++alph_mid) {
					sum_location += newCtc._countValues[alph1][alph_mid][rateCategor] *
						tmp._countValues[alph_mid][alph2][rateCategor];
				}
				resCtc->setCount(alph1,alph2,rateCategor,sum_location);
			}
		}
	}
}



// this function will multiply two CTC's to form a new one.  It assumes that the secount CTC is revesed
void exactSemphyDistance::ctcFromTwoCtcReverse(countTableComponent* resCtc,
											   const countTableComponent& newCtc,
											   const countTableComponent& oldCtc){
  // I'm doing this while assuming that we allways keep the ctc's from the newer nodes to the older ones.
  // this is consistant with computing from node to father, but not the reverse.
  // doing it the other way around whould make more sence computational, as we need to computer the reverse of one 
  // of the ctc's, namely the secound one.  as we allways computer the pairs while using the new CTC, 
  // we can reverse it once and save a lot of work.  not implinted here.
	countTableComponent tmp; 
	tmp.countTableComponentAllocatePlace(_pi.alphabetSize(),_pi.categories());
	tmp.zero();
	//  now, we sum each line to computer p(b)=SIGMA_c (P(b,c)) so we can compute P(c|b)=p(b,c)/p(b)
	int rateCategor =0;
	int alph1 =0;int alph2 =0;
	
	for (rateCategor =0; rateCategor< _pi.categories(); ++rateCategor) {
		for (alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
			MDOUBLE sum(0.0); 
			for (alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
				sum += oldCtc._countValues[alph2][alph1][rateCategor];
			}
			for (alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
				tmp.setCount(alph1, alph2, rateCategor, oldCtc._countValues[alph2][alph1][rateCategor]/(sum+EPSILON));
			}
		}
	}
	for (rateCategor =0; rateCategor< _pi.categories(); ++rateCategor) {
		for (alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
			for (alph2 =0; alph2< _pi.alphabetSize(); ++alph2){
				MDOUBLE sum_location(0.0); 
				for (int alph_mid = 0; alph_mid < _pi.alphabetSize() ; ++alph_mid) { // k
					sum_location += newCtc._countValues[alph1][alph_mid][rateCategor] * tmp._countValues[alph_mid][alph2][rateCategor];
				}
				resCtc->setCount(alph1,alph2,rateCategor,sum_location);
			}
		}
	}
}

void exactSemphyDistance::computeExactDistancesSeparateNodes(
						const tree::nodeP newNode,
						const tree::nodeP oldNode){
	int from=newNode->id();
	int mid=newNode->father()->id();
	int to=oldNode->id();
	if (_tmpvvctcp[from][to]==NULL){	// new 
	    _tmpvvctcp[from][to]=new countTableComponent;
	    _tmpvvctcp[from][to]->countTableComponentAllocatePlace(_pi.alphabetSize(),_pi.categories());
	    _tmpvvctcp[from][to]->zero();
	}
	if (_tmpvvctcp[mid][to]!=NULL){
	  ctcFromTwoCtc(_tmpvvctcp[from][to],*_tmpvvctcp[from][mid],*_tmpvvctcp[mid][to]);
	}
	else {
	  ctcFromTwoCtcReverse((_tmpvvctcp[from][to]),*_tmpvvctcp[from][mid],*_tmpvvctcp[to][mid]);
	}
	acculate_ctc(from,to);
}
*/

