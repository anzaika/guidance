// $Id: minSpanTree.h 409 2005-06-28 13:12:24Z ninio $

#ifndef ___MIN_SPAN_TREE
#define ___MIN_SPAN_TREE

#include "definitions.h"
#include "rearangeTree.h"


using namespace std;

class minSpanTree {
public:
  typedef pair<double,rearrangeTree::intPair> edgePair;
	static rearrangeTree::pairSet span_tree(const VVdouble &weight, double* score);
private:
  // Don't define this function! ??? why
  minSpanTree& operator=(const minSpanTree& rhs);
  //  BNspanTree(const BNspanTree& rhs);
};

#endif

