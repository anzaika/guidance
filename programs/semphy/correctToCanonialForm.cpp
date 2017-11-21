// $Id: correctToCanonialForm.cpp 832 2006-07-25 13:09:38Z ninio $

#include "correctToCanonialForm.h"
#include "logFile.h"
#include <cassert>

correctToCanonialForm::correctToCanonialForm(tree* t1,
											 const VVdouble& distanceTable,
											 const vector<char>& isRealTaxa)
											 :_isRealTaxa(isRealTaxa),_distances(distanceTable) {
	_etPtr = t1;
}

void correctToCanonialForm::getNodeVectorFromId(vector<tree::nodeP>& all) const {
  vector<tree::nodeP> tmp_all;
  _etPtr->getAllNodes(tmp_all,_etPtr->getRoot());
  all.resize(tmp_all.size());
  for(int i=0;i<tmp_all.size();i++)
    all[tmp_all[i]->id()]=tmp_all[i];
}

void correctToCanonialForm::rootToUnrootedTree(tree& et,vector<int> * nodesThatareTakenOut) {
//	CHECK!!!
	LOG(5,<<"0================== rootToUnrootedTree ==================="<<endl); 
	LOGDO(800,et.output(myLog::LogFile(),tree::ANCESTOR));
	LOG(5,<<"0==================================================="<<endl); 
	if (et.getRoot()->getSon(0)->isInternal()) {
		et.getRoot()->getSon(1)->setDisToFather
			( et.getRoot()->getSon(1)->dis2father()+
			et.getRoot()->getSon(0)->dis2father()	);
		et.getRoot()->getSon(0)->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
		// in this case we should fill the 3 sons from the other sons.
		for (int k=0; k < et.getRoot()->getSon(0)->getNumberOfSons(); ++k) {
			et.getRoot()->setSon(et.getRoot()->getSon(0)->getSon(k));
		//	MDOUBLE sumOfDis = et.getRoot()->getSon(1)->getSon(k)->dis2father() + 
		//		et.getRoot()->getSon(1)->dis2father();
		//	et.getRoot()->getSon(1)->getSon(k)->Setdis2father(sumOfDis);
		}
		et.getRoot()->getSon(0)->removeAllSons();
		et.getRoot()->claimSons();
		// there is the chance that the node in son [0] is not a real taxa, and hence should be take out...
//		if (isRealTaxa(et.getRoot()->getSon(0))==false) {
//			nodesThatareTakenOut->push_back(et.getRoot()->getSon(0)->id());
//			et.getRoot()->setSon(et.getRoot()->getSon(et.getRoot()->sons.size()-1),0);
//			et.getRoot()->sons.pop_back();
//		}
	}
	else {// (et.getRoot()->getSon(1)->sons.empty()) 
		et.getRoot()->getSon(1)->setDisToFather
			( et.getRoot()->getSon(0)->dis2father()+
			et.getRoot()->getSon(1)->dis2father()	);
		et.getRoot()->getSon(1)->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
		// in this case we should fill the 3 sons from the other sons.
		for (int k=0; k < et.getRoot()->getSon(1)->getNumberOfSons(); ++k) {
			et.getRoot()->setSon(et.getRoot()->getSon(1)->getSon(k));
			//MDOUBLE sumOfDis = et.getRoot()->getSon(0)->getSon(k)->dis2father() + 
			//	et.getRoot()->getSon(0)->dis2father();
			//et.getRoot()->getSon(0)->getSon(k)->Setdis2father(sumOfDis);
		}
		et.getRoot()->getSon(1)->removeAllSons();
		et.getRoot()->claimSons();
		// there is the chance that the node in son [0] is not a real taxa, and hence should be take out...
//		if (isRealTaxa(et.getRoot()->getSon(1))==false) {
//			nodesThatareTakenOut->push_back(et.getRoot()->getSon(1)->id());
//			et.getRoot()->setSon(et.getRoot()->getSon(et.getRoot()->getNumberOfSons()-1),1);
//			et.getRoot()->sons.pop_back();
//		}

	}

}

void correctToCanonialForm::correctTree() {
	LOG(800,<<"0================== b4 correcting ==================="<<endl); 
	LOGDO(800,_etPtr->output(myLog::LogFile(),tree::ANCESTOR));
	LOG(800,<<"0==================================================="<<endl); 

	vector<int> nodesThatareTakenOut;

	vector<tree::nodeP> all; // this must be called b4 the nodes exclusion...
	getNodeVectorFromId(all);// all[i] is the node with id = i;

	//bool nextStep = false;
	//while (nextStep == false) {
		makeSureRootIsAnHTU();
		recursiveCleanTreeFromExcessNodes(nodesThatareTakenOut,_etPtr->getRoot());
		if (_etPtr->getRoot()->getNumberOfSons()==2) {
			rootToUnrootedTree(*_etPtr,&nodesThatareTakenOut); 
			recursiveCleanTreeFromExcessNodes(nodesThatareTakenOut,_etPtr->getRoot());
		}
		//_etPtr->output(LOG(5,, tree::ANCESTOR,true) );
	//}


	//_etPtr->output(LOG(5,, tree::ANCESTOR,true) );

//	LOG(5,<<"1==================================================="<<endl;//CHEC)K
//	_etPtr->output(LOG(5,, tree::ANCESTOR,true);//CHEC)K
	LOG(800,<<"0================== b4 recursive ==================="<<endl); 
	LOGDO(800,_etPtr->output(myLog::LogFile(),tree::ANCESTOR));
	LOG(800,<<"0==================================================="<<endl); 

	recursiveAddNodesToMisplacedLeaves(nodesThatareTakenOut,_etPtr->getRoot(),all);

//	LOG(5,<<"2==================================================="<<endl;//CHEC)K
//	_etPtr->output(LOG(5,, tree::ANCESTOR,true);//CHEC)K

	recursiveCorrectmOversizedStarsWithExcessNodes(nodesThatareTakenOut,
	  _etPtr->getRoot(),
	  all);  

	deleteAccessNodesAndResetNodeNumbers(nodesThatareTakenOut, all);

//	LOG(5,<<"3==================================================="<<endl;//CHEC)K
//	_etPtr->output(LOG(5,, tree::ANCESTOR,true) ;//CHEC)K
}

//bool correctToCanonialForm::isRealTaxa(const tree::nodeP& n1) {
//	return (_pi->getSeq(n1->id()) != NULL);
//}



//void correctToCanonialForm::makeSureRootHasThreeChildOrMore(){
void correctToCanonialForm::makeSureRootIsAnHTU(){
	vector<tree::nodeP> nodes;
	_etPtr->getAllNodes(nodes,_etPtr->getRoot());
	for (int i=0; i < nodes.size() ;++i) {
		if (_isRealTaxa[nodes[i]->id()] == 0) {
			_etPtr->rootAt(nodes[i]);
			break;
		}
	}
	
	// we now have the root in an HTU
	int numOfRootChild = _etPtr->getRoot()->getNumberOfSons();
	if (numOfRootChild==0) errorMsg::reportError(" error (0) in function makeSureRootIsAnHTU");
	if (numOfRootChild==1) {// out tree looks like n1 has one son x. x has one or more sons...
		if (_isRealTaxa[_etPtr->getRoot()->getSon(0)->id()]) {
			_etPtr->getRoot()->getSon(0)->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
			for (int i=0; i < _etPtr->getRoot()->getSon(0)->getNumberOfSons() ; ++i) {
				_etPtr->getRoot()->setSon(_etPtr->getRoot()->getSon(0)->getSon(i));
			}
			_etPtr->getRoot()->getSon(0)->removeAllSons();
			_etPtr->getRoot()->claimSons();
		}
		else _etPtr->rootAt(_etPtr->getRoot()->getSon(0));
		numOfRootChild = _etPtr->getRoot()->getNumberOfSons();
		assert (numOfRootChild>1);
		//makeSureRootIsAnHTU();
		//return;
	}
	if (numOfRootChild==2) {
		rootToUnrootedTree(*_etPtr);
	}
}

void correctToCanonialForm::recursiveCleanTreeFromExcessNodes(vector<int>& nodesThatareTakenOut, 
									   tree::nodeP inNode){
// this function handle 2 cases:
// 1. there is an HTU that had been positined as a leaf
// 2. there is an HTU that has only one son.
	if (inNode->getNumberOfSons()==0) return;
	for (int i =inNode->getNumberOfSons()-1 ; i >=0 	 ; --i) {
		recursiveCleanTreeFromExcessNodes(nodesThatareTakenOut,inNode->getSon(i));
		if ((inNode->getSon(i)->getNumberOfSons()==0) && (_isRealTaxa[inNode->getSon(i)->id()]==0)) {
			nodesThatareTakenOut.push_back(inNode->getSon(i)->id());
			inNode->setSon(inNode->getSon(inNode->getNumberOfSons()-1),i);
			inNode->removeLastSon();
		} else if ((inNode->getSon(i)->getNumberOfSons()==1)  && !_isRealTaxa[inNode->getSon(i)->id()]) {
			double newDist = inNode->getSon(i)->getSon(0)->dis2father()+inNode->getSon(i)->dis2father();
			nodesThatareTakenOut.push_back(inNode->getSon(i)->id());
			inNode->setSon(inNode->getSon(i)->getSon(0),i);
			inNode->getSon(i)->setFather(inNode);
			inNode->getSon(i)->setDisToFather(newDist);
		}
    }
}


void SwapSon(tree::nodeP fatherP,
			 tree::nodeP oldSon,
			 tree::nodeP newSon) {
	if (fatherP == NULL) return;
	for(int i=0;i<fatherP->getNumberOfSons();i++) {
	
    if (fatherP->getSon(i)==oldSon){
      fatherP->setSon(newSon,i);
      return;
    }
  }
  errorMsg::reportError("The bird in the tree is not sitting on this Branch..");
}



void correctToCanonialForm::recursiveAddNodesToMisplacedLeaves(
									vector<int>& nodesThatareTakenOut,
									tree::nodeP inNode,
									const vector<tree::nodeP>& all) {
//-------------------------------------------------------------------------
// fix nodes that are leaves and are not leaves in the span tree:
// if the leaf has a child, the leaf is replaced with a new node,
// and the leaf is now the child of this new node...
//-------------------------------------------------------------------------

	if (inNode->getNumberOfSons()==0) return;
	for (int i =0 ; i <inNode->getNumberOfSons(); ++i) {
		recursiveAddNodesToMisplacedLeaves(nodesThatareTakenOut,inNode->getSon(i),all);
	}// we start from the leaves towards the root.

	if (_isRealTaxa[inNode->id()]) {// i.e. taxa with childs...
		if (nodesThatareTakenOut.empty()) {
			errorMsg::reportError(" error in function recursiveAddNodesToMisplacedLeaves -  not enough nodes to add. ");
		}
		tree::nodeP fill_node = all[nodesThatareTakenOut.back()];	// get an unused node
		nodesThatareTakenOut.pop_back(); // clear it from the list
	    tree::nodeP old_father=inNode->father();
		fill_node->setFather(old_father);
		fill_node->copySons(inNode);
		fill_node->setSon(inNode);
		fill_node->setDisToFather(inNode->dis2father());
		inNode->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
		inNode->setFather(fill_node);
		inNode->removeAllSons();
	    SwapSon(old_father, inNode, fill_node);
		fill_node->claimSons();
		old_father->claimSons();
		//if (_etPtr->isRoot(inNode)) _etPtr->setAsRootSemphy(inNode->father);
	}
}

void correctToCanonialForm::getBestPair(const vector<int> &sons,
				 int& fromi,
				 int& toi,
				 const bool useNJToBrakeUpStars){
	if (useNJToBrakeUpStars) {
		MDOUBLE min=VERYBIG;
		for (int i=0;i<sons.size()-1;++i){
			for (int j=i+1 ;j<sons.size();++j){
				if (_distances[sons[i]][sons[j]]<min) {
					min=_distances[sons[i]][sons[j]];
					fromi=i;
					toi=j;
				}
			}
		}
	} else {
		fromi=0;
		toi=1;
	}
}


 // this function choses the star brakeup rather arbitrarely, and shold be looked into. 

void correctToCanonialForm::recursiveCorrectmOversizedStarsWithExcessNodes(
								vector<int>& nodesThatareTakenOut,
								tree::nodeP nodeToCorrect,
								const vector<tree::nodeP>& all){
	if (nodeToCorrect->getNumberOfSons()==0) return; // no correction is needed
//	bool node2correctIsTheRoot = _etPtr->isRoot(nodeToCorrect);
	bool node2correctIsTheRoot = nodeToCorrect->isRoot();
	bool correctionIsNeeded = (
			(!node2correctIsTheRoot && 
			nodeToCorrect->getNumberOfSons()>2) || 
			
			(node2correctIsTheRoot && 
			nodeToCorrect->getNumberOfSons()>3)
	);
	if (correctionIsNeeded) {

  		static vector<int> sons;
		sons.clear();
		for (int iSon=0;iSon<nodeToCorrect->getNumberOfSons();iSon++) {
		  sons.push_back(nodeToCorrect->getSon(iSon)->id());
		}
		if (!node2correctIsTheRoot) sons.push_back(nodeToCorrect->father()->id());
		int fromi, toi;
		getBestPair(sons, fromi, toi, true); // true = use NJ to break nodes

		tree::nodeP fill_node = all[nodesThatareTakenOut.back()];
		nodesThatareTakenOut.pop_back();
    
		if (toi==nodeToCorrect->getNumberOfSons()) {	// the father
			fill_node->setFather(nodeToCorrect);
			fill_node->copySons(nodeToCorrect);
			fill_node->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
			nodeToCorrect->removeAllSons();
			nodeToCorrect->setSon(fill_node);
			nodeToCorrect->setSon(fill_node->getSon(fromi));
			fill_node->setSon(fill_node->getLastSon(),fromi); // this will work even if the the fromi is last in the vector
			fill_node->removeLastSon();
			  //      _myBNtree.claimSons(nodeToCorrect->father);
			fill_node->claimSons();
			nodeToCorrect->claimSons();
		} else {			// two diffrent sons.  this time, we join them under us, and then recall this function on our self, and return
			fill_node->setFather(nodeToCorrect);
			fill_node->removeAllSons();
			fill_node->setSon(nodeToCorrect->getSon(fromi));
			fill_node->setSon(nodeToCorrect->getSon(toi));
			fill_node->setDisToFather(VERY_SMALL_DIST_TO_FATHER);
			nodeToCorrect->setSon(fill_node,fromi);
			nodeToCorrect->setSon(nodeToCorrect->getLastSon(),toi); // this will work even if the the toi is last in the vector
			nodeToCorrect->removeLastSon();
			fill_node->claimSons();
			nodeToCorrect->claimSons();
		}
					// now, we may still have too meny children, so we call ourself, and then return      
		recursiveCorrectmOversizedStarsWithExcessNodes(nodesThatareTakenOut,
			nodeToCorrect,all);
		return;
	}  
	for (int i =0 ; i <nodeToCorrect->getNumberOfSons(); ++i) {
		recursiveCorrectmOversizedStarsWithExcessNodes(nodesThatareTakenOut,
		nodeToCorrect->getSon(i),all);
	}
}


void correctToCanonialForm::deleteAccessNodesAndResetNodeNumbers( vector<int>& nodesThatareTakenOut,
																  const vector<tree::nodeP>& all){
  if (nodesThatareTakenOut.size()==0) return;
  int NumNodes= all.size();
  while (nodesThatareTakenOut.size()) {
	tree::nodeP toRemove=all[nodesThatareTakenOut.back()];
	int idOfRemoved = toRemove->id();
	all[NumNodes-1]->setID(idOfRemoved);
	--NumNodes;
	delete toRemove;
	nodesThatareTakenOut.pop_back();
  }
  _etPtr->updateNumberofNodesANDleaves();
}
