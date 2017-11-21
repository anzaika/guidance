// $Id: rearangeTree.cpp 409 2005-06-28 13:12:24Z ninio $

#include "rearangeTree.h"
#include "logFile.h"


rearrangeTree::rearrangeTree(// tree* inT,
					pairSet* inNodeSet,
			     const semphyDistance* inSemphyDis): 
							       _nodesSet(inNodeSet), 
							       _semDist(inSemphyDis)
{};
  
void rearrangeTree::getNodeVectorFromId(vector<tree::nodeP>& all) const{
    vector<tree::nodeP> tmp_all;
	_t1->getAllNodes(tmp_all,_t1->getRoot());
	all.resize(tmp_all.size());
	for(int i=0;i<tmp_all.size();i++) {
		all[tmp_all[i]->id()]=tmp_all[i];
	}
}

void rearrangeTree::reconstructTree(tree& inT) {	// this will change the tree!
  LOG(150,<<"rearrangeTree - enter"<<endl);
  LOGDO(150,inT.output(myLog::LogFile(),tree::ANCESTORID));
  LOG(150,<<"\n\n"<<endl);


  _t1=&inT;
  vector<tree::nodeP> all;
  getNodeVectorFromId(all); // return a vector of nodes such that all[i] is the node with id = i;
  for (vector<tree::nodeP>::iterator i=all.begin(); i != all.end() ;++i) {
    if ((*i)->isLeaf()) {
      _t1->rootAt(*i);
      break;
    }
  }
   recursiveReconstructTree(all, _t1->getRoot(),NULL);

    _t1->rootAt(_t1->getRoot()->getSon(0));	// there has to be a son, and it has to be a non-leaf.
  LOG(150,<<"rearrangeTree - exit"<<endl);
  LOGDO(150,inT.output(myLog::LogFile(),tree::ANCESTORID));
  LOG(150,<<"\n\n"<<endl);

}

void rearrangeTree::recursiveReconstructTree(vector<tree::nodeP>& all,
				       tree::nodeP inNode,
				       tree::nodeP fatherNode) {
//	inNode->_father = fatherNode;
	inNode->setFather(fatherNode);
	inNode->removeAllSons();
	LOG(150, <<"doing "<<inNode->id()<<endl;);	
	pairSet::iterator it = _nodesSet->begin();
	while (it != _nodesSet->end()) {
		// the incrementel part was moved into the function.
	  int nextSonId = -1;
	  bool goodPair = nodeConnectedTo(*it, inNode->id(), nextSonId);
	  LOG(150, <<"check "<<*it<<" got "<<goodPair << endl;);
	  if (goodPair) {
	    inNode->setSon(all[nextSonId]);
	    double tmpDis2father = _semDist->getDistance(inNode->id(),nextSonId);
	    all[nextSonId]->setDisToFather(tmpDis2father);
	    pairSet::iterator tmp(it);
	    ++it;
	    _nodesSet->erase(*tmp); // clear used edges and move forward by one
	    // clear the calues without fucking the itratior up
	  }
	  else ++it;
	}
	//	_nodesSet->remove(intPair(-1,-1)); // clear used edges
	for(int i=0;i<inNode->getNumberOfSons();i++) {
	  recursiveReconstructTree(all,inNode->getSon(i),inNode);
	}
}

bool rearrangeTree::nodeConnectedTo(const intPair& inPair,
				    const int inNodeID, int &nextSonId) {
  //inNodeID is the nodeID we are looking for in the pair
  //next son id is the node it is connected to.
  if (inPair.first == inNodeID) nextSonId = inPair.second;
  else if (inPair.second == inNodeID) nextSonId = inPair.first;
  else return false;
  return true;
}


ostream& operator<< (ostream &sout,  const rearrangeTree::intPair& rpair) {
  sout << "<"<<rpair.first<<","<<rpair.second<<">";
  return sout;
}


