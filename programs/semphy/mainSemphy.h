// $Id: mainSemphy.h 2324 2007-08-08 13:55:45Z privmane $

// This comment is a test for a problem with SVN commit

#ifndef ___MAIN_SEMPHY
#define ___MAIN_SEMPHY

#include "aaJC.h"
#include "alphabet.h"
#include "amino.h"
#include "chebyshevAccelerator.h"
#include "clustalFormat.h"
#include "definitions.h"
#include "distanceMethod.h"
#include "fastaFormat.h"
#include "gammaDistribution.h"
#include "hky.h"
#include "maseFormat.h"
#include "molphyFormat.h"
#include "nucJC.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "readDatMatrix.h"
#include "sequenceContainer.h"
#include "semphy_cmdline.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h"
#include "cmdline2EvolObjs.h"
#include "distanceBasedSeqs2TreeFactory.h"

#include <iostream>
#include <vector>
#include <string>
using namespace std;

class mainSemphy {
public:
	explicit mainSemphy(int argc, char* argv[]);
	explicit mainSemphy(const gengetopt_args_info& gn);


	virtual ~mainSemphy();
	void compute(bool doBootstrapOneIteration=false);
	enum AlgorithmType {DO_DEFAULT = 0, DO_SEMPHY=1, DO_BBL=2, DO_LIKELIHOOD=4, DO_PERPOS_LIKELIHOOD=8, DO_NJ=16};
	tree getTree() const {return *_etPtr;}
	const stochasticProcess& getStochasticProcess() const {return *_spP;}
	MDOUBLE getLikelihood() const { return _treeLogLikelihood;}
	void setWeights(const Vdouble& weights);
	sequenceContainer getSequenceContainer() const { return _sc;}

	void setGlobalRate(MDOUBLE rate) {_spP->setGlobalRate(rate);};
	void output() const;
	ostream& out() const {return *_outPtr;}
	void setTree(const tree& inEtPtr);
	void nullifyTree();
  void semphyCorrectToCanonialForm(void);
	void initializeFromArgsInfo();
	void computeNJtree(bool doBootstrapOneIteration=false);
	void extractPosteriorProb(VVdoubleRep & posteriorProbVV) const;
	void optimizeBranchLengths();
	void computeSemphyTree();
    void optimizeAlphaOnly();  // NOTE: If an SSRV model is used this method also optimizes the Nu parameter
	void optimizeGlobalRate();
	void computeLikelihoodAndLikelihoodPerPosition();
    void computeLikelihoodAndLikelihoodPerPositionAndPosterior();
	void printLikelihoodAndLikelihoodPerPosition() const;
    void printPosterior() const;
	void printTree(const int logLvl=3) const;
	void printTreeToTreeFile() const;
	void computeTree(bool doBootstrapOneIteration=false);
	void optimizeParameters();
	void getDistanceTableAndNames(VVdouble& disTable,
								  vector<string> & vNames,
								  const distanceMethod* cd) const;
	void computeNJtreeFromDisTableAndNames(const VVdouble& disTable,
										  const vector<string> & vNames);

	gengetopt_args_info _args_info; // the information form the command line
    cmdline2EvolObjs<gengetopt_args_info> _evolObj;

private:
	void printSemphyTitle(ostream & out);
	void readCommandLineInformation(int argc, char* argv[]);
	void initializeRandomSeed() const;
	void initializeLogFile() const;
	void initializeAlphabet();
	void takeCareOfGaps();
	void readTreeFile();
	void readConstraintTreeFile();
	void initializeStochaticProcess();
	void initializeSequenceContainer();
	void initializeOutputStream();
	void argsConsistencyCheck() const;
	void constraintTreeConsistencyCheck() const;


	tree* _etPtr;
	tree* _constraintTreePtr;
	ostream* _outPtr;
	Vdouble* _weights;
        distanceBasedSeqs2Tree* _s2tPtr;

	int _numberOfRandomStart;
	VVdoubleRep _posterior; // per position posterior of rates
	VdoubleRep _llpp; // log likelihood per position
	MDOUBLE _treeLogLikelihood; // log likelihood of the tree

	stochasticProcess *_spP;
	sequenceContainer _sc;
	alphabet*  _alphP;
	VVdoubleRep _posteriorRates;
};


#endif


