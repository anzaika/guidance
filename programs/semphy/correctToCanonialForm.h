// 	$Id: correctToCanonialForm.h 922 2006-09-21 09:36:33Z ninio $	

#ifndef ___CORRECT_TO_CANONIAL_FORM
#define ___CORRECT_TO_CANONIAL_FORM

#include "tree.h"
#include "definitions.h"

//const MDOUBLE VERY_SMALL_DIST_TO_FATHER = 0.001;
//const MDOUBLE VERY_SMALL_DIST_TO_FATHER = 1e-9;
const MDOUBLE VERY_SMALL_DIST_TO_FATHER = 0.000001;

class correctToCanonialForm {
public:
	explicit correctToCanonialForm(tree* t1,
					const VVdouble& distanceTable,
					const vector<char>& isRealTaxa);
	void correctTree();
	tree* getEvolTreePtr() const {return _etPtr;}
private:
	void makeSureRootIsAnHTU();
	void rootToUnrootedTree(tree& et,vector<int> * nodesThatareTakenOut = NULL);
//	bool isRealTaxa(const tree::nodeP& n1);
	void recursiveCleanTreeFromExcessNodes(vector<int>& nodesThatareTakenOut, 
									   tree::nodeP inNode);

	void recursiveAddNodesToMisplacedLeaves(
									vector<int>& nodesThatareTakenOut,
									tree::nodeP inNode,
									const vector<tree::nodeP>& all);

	void getBestPair(const vector<int> &sons,
				 int& fromi,
				 int& toi,
				 const bool useNJToBrakeUpStars);

	void recursiveCorrectmOversizedStarsWithExcessNodes(
								vector<int>& nodesThatareTakenOut,
								tree::nodeP nodeToCorrect,
								const vector<tree::nodeP>& all);
	void getNodeVectorFromId(vector<tree::nodeP>& all) const;
	void deleteAccessNodesAndResetNodeNumbers(vector<int>& nodesThatareTakenOut, const vector<tree::nodeP>& all);

private:
	tree *_etPtr;
	const vector<char>& _isRealTaxa; // _isRealTaxa[i] = 0 <==> node with id i is HTU
	const VVdouble& _distances; // distances between any two nodes. consider remove...
//	positionInfo* _pi;  must be here, so that we can know which node is real taxa.
};

#endif

