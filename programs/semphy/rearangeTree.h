// $Id: rearangeTree.h 409 2005-06-28 13:12:24Z ninio $

#ifndef ___REARANGE___TREE
#define ___REARANGE___TREE

#include <set>
using namespace std;

#include "semphyDistance.h"
#include "tree.h"


class rearrangeTree {
public:
  	typedef pair<int,int> intPair;
        typedef set<intPair> pairSet;

	rearrangeTree(
						pairSet* inNodeSet,
						const semphyDistance* aSD);

	void reconstructTree(tree& t1);
private:

	void recursiveReconstructTree(vector<tree::nodeP>& all,
						   tree::nodeP inNode,
						   tree::nodeP fatherNode);

	bool nodeConnectedTo(const intPair& inPair,
					  const int inNodeID, int &nextSonId);

	void getNodeVectorFromId(vector<tree::nodeP>& all) const;
	
	
private:
	tree* _t1;
	pairSet* _nodesSet;
	const semphyDistance* _semDist;
};

ostream& operator<< (ostream &sout,  const rearrangeTree::intPair& rpair);


#endif

