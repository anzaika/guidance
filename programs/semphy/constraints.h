
// 	$Id: constraints.h 922 2006-09-21 09:36:33Z ninio $	
// general Idea:  

// 1. you start off by reading a constrain- basically a tree.

// 2. you change the root of the original tree to fit the root of the
//    constraint tree.

// 3. you run over the original tree and over the constraint tree, and for
//    each internal node you list the names of all the leafs under that
//    node.  for the original tree you collect the internal nodes as well
//    (including the node itself)

// 4. you go over each of the constraints, chose a leaf from the set, go to
//    that leaf and clime up the tree until you find a node with exactly
//    the same set of leafs.  you then activate the constraint for the set
//    of all nodes under that node (from the list)

// 5. return the constraint matrix.

// NOTE: the same constraint matrix may need to be recomputed after each
// semphy iteration - the numbers of the internal nodes may be shifted
// around during the semphy step.


#ifndef ___CONSTRAINTS
#define ___CONSTRAINTS

#include <string>
#include <set>
#include <vector>
#include "tree.h"

using namespace std;

class constraints {
public: 
  explicit constraints(const string& constraingFileName);
  explicit constraints(const tree& in_ct);

  void setTree(const tree& t1);	// give the current phlogony, recompute penerly table
  const VVdouble& getPeneltyTable();
  bool fitsConstraints();
  void outputMissingClads(ostream& out);
  
  
  typedef set<string> childrenSet;
  typedef set<childrenSet> childrenSetSet;
  typedef pair<childrenSet, set<int> > csidp;

private:
  tree _ct;			// constraints
  tree _currentTree;
  VVdouble _peneltyTable;
  string _sonOfRoot;
  childrenSetSet _css, _cssTmp;
  
  void constraints_real_constractor();
  csidp recursivlyFillTable(const tree::nodeP root);
  void addPenetly(set<int>& v);
  void make_list_of_constraints();
  childrenSet recursive_make_list_of_constraints(const tree::nodeP & root);
};

#endif // ___CONSTRAINTS
