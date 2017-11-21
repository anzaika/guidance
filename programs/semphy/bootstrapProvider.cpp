// $Id: bootstrapProvider.cpp 1282 2006-12-05 11:03:10Z privmane $

#include "bootstrapProvider.h"
#include "getRandomWeights.h"
#include "logFile.h"
#include <ostream>
using namespace std;

bootstrapProvider::bootstrapProvider(const gengetopt_args_info& gn) :_bp(NULL) {
	_args_info = gn;
}

void bootstrapProvider::createRandomWeights(const int seqLen){
	_weights.resize(seqLen,0.0);
	getRandomWeights::standardBPWeights(_weights);
}

void bootstrapProvider::output(ostream & out) const {
	if (_args_info.BPconsensus_given) {
		out<<"# The consensus tree is: "<<endl;
		_consensusTree.output(out);

		LOG(3,<<"# The consensus tree is: "<<endl);
		LOGDO(1,_consensusTree.output(out));
	}
	//2. print the weights on the user tree.
		out<<"# The bootstrap values imposed on the resulting tree: "<<endl;
		LOG(3,<<"# The bootstrap values imposed on the resulting tree: "<<endl);
		_bp->printTreeWithBPvalues(out,_inputTree,_treeSupport);
		_bp->printTreeWithBPvalues(myLog::LogFile(),_inputTree,_treeSupport);
}


void bootstrapProvider::computeBP(mainSemphy & ms) {
	_inputTree = ms.getTree();
	const int numberOfBP = _args_info.BPrepeats_arg;
	int seqLen = ms.getSequenceContainer().seqLen();
	for (int i=0; i < numberOfBP; ++i) {
		createRandomWeights(seqLen);
		ms.nullifyTree();
		ms.setWeights(_weights);
		ms.compute(true);
		_treeVec.push_back(ms.getTree());
		_likelihoodVec.push_back(ms.getLikelihood());
		_stochasticProcessVec.push_back(ms.getStochasticProcess());
	}
	if (_bp) delete _bp;
	_bp = new bootstrap(_treeVec);

	if (_args_info.BPconsensus_given) {
		computeConsensus();
	} 
    computeTreeSupport(_inputTree);
}

void bootstrapProvider::computeTreeSupport(const tree& et) {
	
	_treeSupport = _bp->getWeightsForTree(et);
}

void bootstrapProvider::computeConsensus() {
	const int treshold = _args_info.BPconsensus_arg;
	computeConsensus(treshold);
}
	

void bootstrapProvider::computeConsensus(const MDOUBLE treshold) {
	_consensusTree = _bp->consensusTree();
}
