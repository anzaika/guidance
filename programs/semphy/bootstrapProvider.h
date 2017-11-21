// 	$Id: bootstrapProvider.h 1282 2006-12-05 11:03:10Z privmane $	

#ifndef ___BOOTSTRAP_PROVIDER
#define ___BOOTSTRAP_PROVIDER

#include "bootstrap.h"
#include "definitions.h"
#include "mainSemphy.h"
#include "semphy_cmdline.h"
#include "stochasticProcess.h"
#include "tree.h"
#include <vector>
#include <map>
using namespace std;

class bootstrapProvider {
public:
	explicit bootstrapProvider(const gengetopt_args_info& gn);
	
	// this function creates the many tree.
	void computeBP(mainSemphy & ms);
	
	void computeConsensus(const MDOUBLE treshold);
	void computeConsensus();
	void output(ostream & out) const;


private:
	void createRandomWeights(const int seqLen);
	void computeTreeSupport(const tree& et);

	gengetopt_args_info _args_info;
	vector<tree> _treeVec;
	vector<MDOUBLE>	_likelihoodVec;
	vector<stochasticProcess> _stochasticProcessVec;
	vector<MDOUBLE> _weights;
	tree _consensusTree;
	bootstrap* _bp;

	tree _inputTree; // the tree on which support is computed.
	map<int, MDOUBLE> _treeSupport;
};

#endif


