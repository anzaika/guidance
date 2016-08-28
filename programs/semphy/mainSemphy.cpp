// $Id: mainSemphy.cpp 6002 2009-03-20 19:39:03Z privmane $

// DO: (1) ADD NUMBER OF DESCRETE CATEGORIES TO THE OPTIONS.
// DO: (2) ADD NUMBER OF RANDOM STARTS TO THE OPTIONS.
// DO: (3) CHENYSHEV PARAMETERS TO THE OPTIONS.
// DO: (4) ADD PARAMETER CONCERNING THE STARTING NJ TREE.
// DO: (5) ADD PARAMETER CONCERNING THE GAMMA PARAM OPTIMIZATION.

#include "definitions.h"
#include "mainSemphy.h"
#include "jcDistance.h" 
#include "distanceTable.h" 
#include "nj.h" 
#include "constraints.h"
#include "bblEM.h"
#include "semphySearchBestTree.h"
#include "logFile.h"
#include "readDatMatrix.h"
#include "datMatrixHolder.h"
#include "likelihoodComputation.h"
#include "likeDist.h"
#include "bestAlpha.h"
#include "bestParamUSSRV.h"
#include "bestTamura92param.h"
#include "bestGtrModelParams.h"
#include "findRateOfGene.h"
#include "talRandom.h"
#include "someUtil.h"
#include "getRandomWeights.h"
#include "codon.h"
#include "recognizeFormat.h"
#include "generalGammaDistributionLaguerre.h"
#include "correctToCanonialForm.h"
#include "ssrvDistanceSeqs2Tree.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdio>
using namespace std;

//*******************************************************************************
// constructors
//*******************************************************************************
mainSemphy::mainSemphy(const gengetopt_args_info& gn):	_args_info(gn),_evolObj(_args_info),_weights(NULL), _s2tPtr(NULL), _numberOfRandomStart(1){
	initializeFromArgsInfo();
}

mainSemphy::mainSemphy(int argc, char* argv[]):_weights(NULL), _s2tPtr(NULL) , _numberOfRandomStart(1){
	readCommandLineInformation(argc,argv);
	initializeFromArgsInfo();
	myLog::printArgv(1, argc, argv);

}

void mainSemphy::readCommandLineInformation(int argc, char* argv[]) {
	if (argc < 2) errorMsg::reportError("The program must get some parameters in the command line, use -h for help");
	if (cmdline_parser(argc, argv, &_args_info) != 0) {
	  errorMsg::reportError("error reading command line",1);
	}	
	cmdline2EvolObjs<gengetopt_args_info> _evolObj(_args_info);
}

void mainSemphy::initializeFromArgsInfo(){
//	AND WE CHECK TO SEE IF WE ARE ASKED FOR SOMETING WE CAN NOT DELIVER
//	if	(_args_info.min_improv_given)  errorMsg::reportError("minimum-improvement not yet implimented");
//	if	(_args_info.maxDistance_given) errorMsg::reportError("max distance not yet implimented");
//	if	(_args_info.exact_given)		  errorMsg::reportError("exact counts not yet implimented");
	argsConsistencyCheck();
	_evolObj.initializeRandomSeed();
	_evolObj.initializeLogFile();
	//   initializeAlphabet();
	_alphP=_evolObj.cmdline2Alphabet();
	//	initializeSequenceContainer();
	_sc = _evolObj.cmdline2SequenceContainer(_alphP);
	//	takeCareOfGaps();
	_evolObj.takeCareOfGaps(_sc);
	//	readTreeFile();
	_etPtr=_evolObj.cmdline2Tree();
	//  make sure tree is in canonical form
    semphyCorrectToCanonialForm();
	//	readConstraintTreeFile();
	_constraintTreePtr = _evolObj.cmdline2ConstraintTree();
	//	initializeStochaticProcess();
	_spP=_evolObj.cmdline2StochasticProcessSafe();
	//	initializeOutputStream(); //out().
	_outPtr = _evolObj.cmdline2OutputStream();
	printSemphyTitle(out());
	constraintTreeConsistencyCheck(); // check that the tree is compatible with the constraint tree if given.
	if (_args_info.posteriorRates_given)
		_posteriorRates = _evolObj.cmdline2PosteriorRates();
}

void mainSemphy::printSemphyTitle(ostream & out) {
  if (_args_info.verbose_arg<0) return;	// just don't print if using -v-1
	out<<"#################################################################"<<endl;
	out<<"# SEMPHY: A Structural EM Algorithm for Phylogenetic Inference  #"<<endl;
	out<<"# for information, please send email to semphy@cs.huji.ac.il    #"<<endl;
	out<<"#################################################################"<<endl;
	out<<endl;
}

mainSemphy::~mainSemphy() {
	if (_alphP) delete _alphP;
	if (_spP) delete _spP;
	if (_etPtr) delete _etPtr;
	if (_constraintTreePtr) delete _constraintTreePtr;
	if (_weights) delete _weights;
	if (_s2tPtr) delete _s2tPtr;
	myLog::endLog();// close log file as late as posible
}
void mainSemphy::semphyCorrectToCanonialForm(void){
  if (_etPtr == NULL) return;

  int nNodes=_etPtr->getNodesNum();

  vector<char> isRealTaxa(nNodes,0);
  vector<tree::nodeP> all;
  _etPtr->getAllNodes(all,_etPtr->getRoot());
  for (vector<tree::nodeP>::iterator i=all.begin();i!=all.end();++i)
	isRealTaxa[(*i)->id()]=(*i)->isLeaf();

  VVdouble dummyDistanceTable(nNodes);
  for (int i=0;i<nNodes;++i) {
	dummyDistanceTable[i].resize(nNodes);
	mult(dummyDistanceTable[i], 0.0); 	// all values are 0 now;
  }
  correctToCanonialForm ctcf(_etPtr, dummyDistanceTable, isRealTaxa);
  ctcf.correctTree();
}

void mainSemphy::argsConsistencyCheck() const 
{
	// make sure that if the user askes for an action that requires a
	// tree, either he gives a tree or he askes SEMPHY to make one
	if (!(_args_info.tree_given || _args_info.DistanceTableEstimationMethod_group_counter!=0 || _args_info.SEMPHY_given))	// we neither have a tree of can make one
		errorMsg::reportError("Must ask to build a tree by some method or give a tree as input");
	
	// Check that tree is not given for non-iterative NJ
	if (_args_info.tree_given && (_args_info.NJ_given || _args_info.homogeneousRatesDTME_given || _args_info.pairwiseGammaDTME_given))
		errorMsg::reportError("Don't give initial tree with non-iterative NJ");
	
	// Check options related to a given tree and to bootstrap
	if (_args_info.BPonUserTree_given && !_args_info.tree_given)
		errorMsg::reportError("Must give --tree with --BPonUserTree");
    if (_args_info.BPonUserTree_given && _args_info.commonAlphaDTME_given
		&& !_args_info.alpha_given)
		errorMsg::reportError("Must give --alpha with --BPonUserTree when using --commonAlphaDTME",1);
    if (_args_info.BPonUserTree_given && _args_info.rate4siteDTME_given)
		errorMsg::reportError("No support for --BPonUserTree when using --rate4siteDTME",1);
    if (_args_info.BPonUserTree_given && _args_info.posteriorDTME_given
		&& !_args_info.posteriorRates_given)
		errorMsg::reportError("Must give --posteriorRates with --BPonUserTree when using --posteriorDTME",1);
    if (_args_info.BPonUserTree_given && _args_info.optimizeAlpha_given && _args_info.posteriorDTME_given)
		errorMsg::reportError("Don't ask for --optimizeAlpha when using --BPonUserTree with iterative NJ methods - the alpha or rates should be given as input by --alpha or --posteriorRates",1);

	// Check that the method chosen is compatible with the configuration of the other parameters and flags
	if (_args_info.commonAlphaDTME_given || _args_info.posteriorDTME_given) {
		if (!_args_info.optimizeAlpha_given && !_args_info.BPonUserTree_given)
			errorMsg::reportError("Must use --optimizeAlpha with --commonAlphaDTME or --posteriorDTME",1);
	}
	if (_args_info.posteriorRates_given && ! _args_info.posteriorDTME_given)
		errorMsg::reportError("--posteriorRates can only be used with --posteriorDTME");


	if (_args_info.homogeneous_flag && (_args_info.optimizeAlpha_given || _args_info.alpha_given || _args_info.categories_given || _args_info.laguerre_flag)) // problematic with homogeneous
		errorMsg::reportError("Must not use --optimizeAlpha or --alpha with --homogeneous",1);

	if (_args_info.rate4siteDTME_given || _args_info.homogeneousRatesDTME_given) {
		if (!_args_info.homogeneous_flag) {
			errorMsg::reportError("Must use --homogeneous with --rate4siteDTME or --homogeneousRatesDTME_given");
		}
	}

	// Check consistency of options with using the SSRV model
	if (_args_info.ssrv_flag) {
		if (_args_info.homogeneous_flag)
			errorMsg::reportError("Cannot use SSRV with --homogeneous.  Only a Gamma-ASRV model is allowed");
		if (!(_args_info.commonAlphaDTME_given && _args_info.optimizeAlpha_given))
			errorMsg::reportError("Currently, SSRV is only implemented to run with --commonAlphaDTME and --optimizeAlpha");
		if (_args_info.SEMPHY_flag)
			errorMsg::reportError("Currently, SSRV is only implemented to run with --commonAlphaDTME and --optimizeAlpha, and not with --SEMPHY");
	}
}

// check that the inputed tree, constraint tree
// are consistant, to prevent surprises and horrible problems.
void mainSemphy::constraintTreeConsistencyCheck() const {
	if (_etPtr!=NULL) {
		if (_constraintTreePtr !=NULL) {
			constraints c1(*_constraintTreePtr);
			c1.setTree(*_etPtr);
			if (!c1.fitsConstraints()) {
				LOG(1,<<"Input tree does not fit constraints!"<<endl);
			    LOGDO(1,c1.outputMissingClads(myLog::LogFile()));
				errorMsg::reportError("Please enter a starting tree that fits the constraints");
			}
		}
	}
}

void mainSemphy::optimizeBranchLengths() {
	if (_etPtr == NULL) {
		errorMsg::reportError("mainSemphy::optimizeBranchLengths: A tree must be given before optimizing branch length");
	}
	if(_args_info.optimizeAlpha_flag) {
		MDOUBLE bestAlpha = -1;
		if (dynamic_cast<tamura92*>(_spP->getPijAccelerator()->getReplacementModel())) {
			// optimizing params of the tamura92 model
			bestTamura92ParamAlphaAndBBL tmpbestAlpha(*_etPtr,_sc,*_spP, _weights, 5, 0.05,
													  _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, 
													  _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, 
													  _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, 
													  _args_info.epsilonLikelihoodImprovement4BBL_arg, 
													  5.0, _args_info.maxNumOfBBLIter_arg, 1.5 );
			bestAlpha = tmpbestAlpha.getBestAlpha();
			_treeLogLikelihood = tmpbestAlpha.getBestL();


		} else if (dynamic_cast<gtrModel*>(_spP->getPijAccelerator()->getReplacementModel())) {
			// optimizing params of the gtr model
			bestGtrModel optimizer(*_etPtr, _sc, *_spP, _weights, 5,
								   _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg,
								   _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg,
								   true, true);
			bestAlpha=optimizer.getBestAlpha();
			_treeLogLikelihood = optimizer.getBestL();

		} else {
			bestAlphaAndBBL tmpbestAlpha(*_etPtr,_sc,*_spP, _weights, 1.5, 5.0,
										 _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, 
										 _args_info.epsilonLikelihoodImprovement4BBL_arg, 
										 _args_info.maxNumOfBBLIter_arg);
			bestAlpha = tmpbestAlpha.getBestAlpha();
			_treeLogLikelihood = tmpbestAlpha.getBestL();
		}
		
		out()<<"# Best alpha after branch length optimiziation"<<endl;
		out()<<bestAlpha <<endl;
		out()<<"# The likelihood of the tree"<<endl;
		out() <<_treeLogLikelihood<<endl;

		LOG(1,<<"# Best alpha after branch length optimiziation"<<endl);
		LOG(1,<<bestAlpha <<endl);
		LOG(1,<<"# The likelihood of the tree"<<endl);
		LOG(1,<<_treeLogLikelihood<<endl);
	} else {
	  bblEM bblEM1(*_etPtr,_sc,*_spP,_weights,_args_info.maxNumOfBBLIter_arg,_args_info.epsilonLikelihoodImprovement4BBL_arg);//maxIterations=1000,epsilon=0.05
		_treeLogLikelihood=bblEM1.getTreeLikelihood();
		out()<<"# The likelihood of the tree"<<endl;
		out()<<_treeLogLikelihood<<endl;
		LOG(1,<<"# The likelihood of the tree"<<endl);
		LOG(1,<<_treeLogLikelihood<<endl);
	}
}

void mainSemphy::printTree(const int logLvl) const {
	int printMsg=max(3,logLvl);
	out()<<"# The tree"<<endl;
	_etPtr->output(out());
	LOG(printMsg,<<"# The tree"<<endl);
	LOGDO(logLvl,_etPtr->output(myLog::LogFile()));
}

void mainSemphy::printTreeToTreeFile() const {
	if (_args_info.treeoutputfile_given) {
		ofstream treeO(_args_info.treeoutputfile_arg);
		if (! treeO.is_open()) {
			errorMsg::reportError("can not open tree output file");
		}
		_etPtr->output(treeO);
		treeO.close();
	}
}

void mainSemphy::nullifyTree() {
	if (_etPtr) delete (_etPtr);
	_etPtr = NULL;
}

// This function is used, so that for example in bp, a new tree will be computed in NJ.
// The computeNJ for example will not compute the NJ tree if a tree is given.
void mainSemphy::setTree(const tree& inEtPtr) {
	if (_etPtr) delete (_etPtr);
	_etPtr = new tree(inEtPtr);
}

// void mainSemphy::getDistanceTableAndNames(VVdouble& disTable,
// 										  vector<string> & vNames,
// 										  const distanceMethod* cd) const {

// 	giveDistanceTable(cd,_sc,disTable,vNames,_weights);
// }

// void mainSemphy::computeNJtreeFromDisTableAndNames(const VVdouble& disTable,
// 										  const vector<string> & vNames) {
// 	NJalg nj1;
// 	if (_args_info.constraint_given) { // did we get a constraint tree
// 		setTree(nj1.computeTree(disTable,vNames,_constraintTreePtr));
// 	} else {
// 		setTree(nj1.computeTree(disTable,vNames));
// 	}
// }


void mainSemphy::computeNJtree(bool doBootstrapOneIteration) {

	distanceBasedMethod_t dtme = homogeneousRatesDTME; // default, vanila NJ
	if (_args_info.homogeneousRatesDTME_given || _args_info.NJ_given)
		dtme = homogeneousRatesDTME;
	else if (_args_info.pairwiseGammaDTME_given)
		dtme = pairwiseGammaDTME;
	else if (_args_info.commonAlphaDTME_given)
		dtme = commonAlphaDTME;
	else if (_args_info.rate4siteDTME_given)
		dtme = rate4siteDTME;
	else if (_args_info.posteriorDTME_given)
		dtme = posteriorDTME;
	else if (_args_info.SEMPHY_given)
		dtme = homogeneousRatesDTME;
	else errorMsg::reportError("mainSemphy::computeNJtree: An unsuppored DTME was specified");

	bool useJcDistance = (_args_info.nucjc_given || _args_info.aaJC_given);

	if (!_s2tPtr)
		_s2tPtr = distanceBasedSeqs2TreeFactory(dtme, *_spP, useJcDistance, _args_info.optimizeAlpha_flag, _args_info.ssrv_flag, _args_info.epsilonLikelihoodImprovement4iterNJ_arg, _args_info.epsilonLikelihoodImprovement4pairwiseDistance_arg, _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, _args_info.epsilonLikelihoodImprovement4BBL_arg, _args_info.maxNumOfBBLIter_arg);

	if (!doBootstrapOneIteration) {

		// No given initial tree
		if (_etPtr == NULL) {
			if (_args_info.commonAlphaDTME_given) {
				if (!_args_info.ssrv_flag) {
					commonAlphaDistanceSeqs2Tree *caS2tPtr = static_cast<commonAlphaDistanceSeqs2Tree*>(_s2tPtr);
					if (_args_info.alpha_given) { // use the given alpha
						setTree(caS2tPtr->seqs2TreeIterative(_sc, _args_info.alpha_arg, _weights, _constraintTreePtr));
					} else {		       // homogeneous rates in first iteration
						setTree(caS2tPtr->seqs2TreeIterative(_sc, _weights, _constraintTreePtr));
					}
				} else { // Using an SSRV model - run with alpha & nu parameters
					ssrvDistanceSeqs2Tree *ssrvS2tPtr = static_cast<ssrvDistanceSeqs2Tree*>(_s2tPtr);
					if (_args_info.alpha_given) { // use the given alpha & nu
						setTree(ssrvS2tPtr->seqs2TreeIterative(_sc, _args_info.alpha_arg, _args_info.nu_arg, _weights, _constraintTreePtr));
					} else {		       // homogeneous rates in first iteration
						setTree(ssrvS2tPtr->seqs2TreeIterative(_sc, _weights, _constraintTreePtr));
					}
				}

			} else if (_args_info.posteriorRates_given) { // posteriorDTME with given initial posteriorRates (from input file)
				posteriorDistanceSeqs2Tree *posteriorS2tPtr = static_cast<posteriorDistanceSeqs2Tree*>(_s2tPtr);
				setTree(posteriorS2tPtr->seqs2TreeIterative(_sc, _args_info.alpha_arg, _posteriorRates, _weights, _constraintTreePtr));

			} else {  // all other methods
				setTree(_s2tPtr->seqs2Tree(_sc, _weights, _constraintTreePtr));
			}

		// An initial tree (--tree) was given so pass it to the iterative seqs2Tree method
		// NOTE: argsConsistencyCheck makes sure that non-interative NJ can't be run with --tree
		} else {
			if (!_args_info.ssrv_flag) {
				iterativeDistanceSeqs2Tree *itS2tPtr = static_cast<iterativeDistanceSeqs2Tree*>(_s2tPtr);
				if (_args_info.alpha_given) { // use the given alpha
					if (! _args_info.posteriorRates_given) {
						setTree(itS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _args_info.alpha_arg, _weights, _constraintTreePtr));
					} else { // posteriorDTME with given initial posteriorRates (from input file)
						posteriorDistanceSeqs2Tree *posteriorS2tPtr = static_cast<posteriorDistanceSeqs2Tree*>(_s2tPtr);
						setTree(posteriorS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _args_info.alpha_arg, _posteriorRates, _weights, _constraintTreePtr));
					}
				} else {
					setTree(itS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _weights, _constraintTreePtr));
				}
			} else { // Using an SSRV model - run with alpha & nu parameters
				ssrvDistanceSeqs2Tree *ssrvS2tPtr = static_cast<ssrvDistanceSeqs2Tree*>(_s2tPtr);
				if (_args_info.alpha_given) { // use the given alpha & nu
					if (_args_info.nu_given)
						setTree(ssrvS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _args_info.alpha_arg, _args_info.nu_arg, _weights, _constraintTreePtr));
					else
						setTree(ssrvS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _args_info.alpha_arg, _weights, _constraintTreePtr));
				} else {
					setTree(ssrvS2tPtr->seqs2TreeIterative(_sc, *_etPtr, _weights, _constraintTreePtr));
				}
			}
		}

	// Do one bootstrap iteration

	// If initial alpha or nu (for using an SSRV model) were given as
	// commandline input then they were already given to the _s2tPtr (or
	// internal objects) in its construction
	} else {
		if (!_args_info.BPonUserTree_given) {
			// Running bootstrap for the tree constructed during this run
			setTree(_s2tPtr->seqs2TreeBootstrap(_sc, _weights, _constraintTreePtr));
		} else {
			// Running bootstrap on a given user tree:
			// If we use an iterative distance method (commonAlpha or posterior)
			// then side info must have been given as input too
			// so we need to pass it as argument
			if (dtme == commonAlphaDTME) {
				commonAlphaDistanceSeqs2Tree *caS2tPtr = static_cast<commonAlphaDistanceSeqs2Tree*>(_s2tPtr);
				setTree(caS2tPtr->seqs2TreeBootstrap(_sc, _args_info.alpha_arg, _weights, _constraintTreePtr));
			} else if (posteriorDTME) {
				posteriorDistanceSeqs2Tree *posteriorS2tPtr = static_cast<posteriorDistanceSeqs2Tree*>(_s2tPtr);
				setTree(posteriorS2tPtr->seqs2TreeBootstrap(_sc, _posteriorRates, _weights, _constraintTreePtr));
			} else {
				setTree(_s2tPtr->seqs2TreeBootstrap(_sc, _weights, _constraintTreePtr));
			}
		}
	}
}

void mainSemphy::computeSemphyTree() {
	if (_etPtr == NULL) computeNJtree();
	semphySearchBestTree(_sc,*_etPtr,_constraintTreePtr,*_spP,out(),_numberOfRandomStart,
						 _args_info.optimizeAlpha_flag, 
						 _args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg, 
						 _args_info.epsilonLikelihoodImprovement4BBL_arg, 
						 _args_info.maxNumOfBBLIter_arg);
}

void mainSemphy::optimizeGlobalRate() {
    out()<<"we are in void mainSemphy::optimizeGlobalRate()"<<endl;
	// THIS FUNCTION SHOULD BE FIXED IN THE FUTURE TO DO ITERATIONS OVER THE TWO COMPUTATIONS...
	if (_args_info.optimizeAlpha_flag) { // here we do alpha but NO bbl.
		optimizeAlphaOnly();
	}
	MDOUBLE rateOfGene=findTheBestFactorFor(*_etPtr,
											_sc,
											*_spP, // changes *_spP
											_weights,
											_treeLogLikelihood);
	out()<<"# The global rate of the tree"<<endl;
	out() << rateOfGene<<endl;
	out()<<"# The likelihood of the tree is "<<_treeLogLikelihood<<endl;
	LOG(1,<<"# The global rate of the tree"<<endl);
	LOG(1,<< rateOfGene<<endl);
	LOG(1,<<"# The likelihood of the tree is "<<_treeLogLikelihood<<endl);
}

void mainSemphy::optimizeAlphaOnly() {
	if (!_args_info.ssrv_flag) {
	    bestAlphaFixedTree bestAlpha(*_etPtr, _sc, *_spP, _weights, 1.5);
		out()<<"# Best alpha (for fixed branch lengths)"<<endl;
		out()<<bestAlpha.getBestAlpha() <<endl;
		LOG(1,<<"# Best alpha (for fixed branch lengths)"<<endl);
		LOG(1,<<bestAlpha.getBestAlpha() <<endl);
	    _treeLogLikelihood = bestAlpha.getBestL();

	} else {
		// Using SSRV optimizations
		bestParamSSRV* optimizer = NULL;
		if (_args_info.tamura92_given)
			optimizer = new bestParamSSRV(true,true,false,false); // optimize alpha and nu, not tamura92, not bbl		
		else 
			optimizer = new bestParamSSRV(true,true,true,false); // optimize alpha and nu, tamura92, not bbl
		(*optimizer)(*_etPtr,_sc,*(static_cast<stochasticProcessSSRV*>(_spP)),_weights,
					 15,15,0.5,_args_info.epsilonLikelihoodImprovement4alphaOptimiz_arg,_args_info.epsilonLikelihoodImprovement4iterNJ_arg,
					 _args_info.epsilonLikelihoodImprovement4BBL_arg, _args_info.maxNumOfBBLIter_arg);
		out()<<"# Best alpha (for fixed branch lengths)"<<endl;
		out()<<optimizer->getBestAlpha() <<endl;
		LOG(1,<<"# Best alpha (for fixed branch lengths)"<<endl);
		LOG(1,<<optimizer->getBestAlpha() <<endl);
		out()<<"# Best Nu (for fixed branch lengths)"<<endl;
		out()<<optimizer->getBestNu() <<endl;
		LOG(1,<<"# Best Nu (for fixed branch lengths)"<<endl);
		LOG(1,<<optimizer->getBestNu() <<endl);
		_treeLogLikelihood = optimizer->getBestL();
		delete optimizer;
	}

	out()<<"# The likelihood of the tree"<<endl;
	out() <<_treeLogLikelihood<<endl;
	LOG(1,<<"# The likelihood of the tree"<<endl);
	LOG(1,<<_treeLogLikelihood<<endl);
}

void mainSemphy::computeLikelihoodAndLikelihoodPerPosition() {
	_treeLogLikelihood = 0.0;
	_llpp.clear();
	computePijGam cpij;
	cpij.fillPij(*_etPtr,*_spP);
	for (int pos=0; pos < _sc.seqLen() ;++pos) {
		MDOUBLE tmpLL = log(likelihoodComputation::getLofPos(pos,*_etPtr,_sc,cpij,*_spP));
		_treeLogLikelihood += tmpLL;
		_llpp.push_back(tmpLL);
	}
}


void mainSemphy::computeLikelihoodAndLikelihoodPerPositionAndPosterior() {
	_treeLogLikelihood = 0.0;
	_llpp.clear();
	_posterior.clear();
	VdoubleRep posPost(_spP->categories());
	//	getPosteriorOfRatesAndLLPP(*_etPtr, _sc, *_spP, cup, 
	computePijGam cpij;
	cpij.fillPij(*_etPtr,*_spP);
	for (int pos=0; pos < _sc.seqLen() ;++pos) {
	  MDOUBLE tmpLL = log(likelihoodComputation::getLofPosAndPosteriorOfRates(pos,*_etPtr,_sc,cpij,*_spP, posPost));
		_treeLogLikelihood += tmpLL;
		_llpp.push_back(tmpLL);
		_posterior.push_back(posPost);
	}
}

void mainSemphy::printLikelihoodAndLikelihoodPerPosition() const {
	out()<<"# The log likelihood of the tree is:"<<endl;
	out()<<_treeLogLikelihood<<endl;
	LOG(3,<<"# The log likelihood of the tree is:"<<endl);
	LOG(1,<<_treeLogLikelihood<<endl);

	if (_args_info.PerPosLike_given){
		out()<<"# The log likelihood per position:"<<endl;
		out()<<_llpp<<endl;
		LOG(3,<<"# The log likelihood per position:"<<endl);
		LOG(1, <<_llpp<<endl);
	} 
}

void  mainSemphy::printPosterior() const {
	if (_args_info.PerPosPosterior_given){  
	  out()<<"# The posterior of the rates is:"<<endl;
	  out()<<_posterior<<endl;
	  LOG(8,<<"# centroieds for the rate"<<endl);
	  for (int i=0;i<_spP->categories();++i)
		LOG(8,<<"  "<< _spP->rates(i));
	  LOG(8,<< endl);
	  LOG(3,<<"# The posterior of the rates is:"<<endl);
	  LOG(1,<<_posterior<<endl);
	}
}

void mainSemphy::setWeights(const Vdouble& weights) {
	if (_weights) delete _weights;
	_weights = new Vdouble(weights);
}


void mainSemphy::computeTree(bool doBootstrapOneIteration) {
	if (_args_info.SEMPHY_given) {
		computeSemphyTree(); // Note that computeSemphyTree calls computeNJtree
	} else if (_args_info.DistanceTableEstimationMethod_group_counter) { // Do NJ of some sort without SEMPHY
		computeNJtree(doBootstrapOneIteration);
	}
}
 
		
	
void mainSemphy::optimizeParameters() {
	// Don't optimize parameters in case SEMPHY was not called, and we
	// ran an iterative NJ method that already did the optimization
	if (!_args_info.SEMPHY_given
		&& (_args_info.commonAlphaDTME_given || _args_info.posteriorDTME_given))
		return;

	// 1. OPTIMIZE BRANCH LENGTH
	// 2. OPTIMIZE ALPHA
	// 3. OPTIMIZE BOTH
	if (_args_info.bbl_given) {// this will optimize both branch lengths and alpha if needed.
		optimizeBranchLengths();
	} else if (_args_info.rate_flag) { // DO GLOBAL RATE OPTIMIZATION (WITH OR WITHOUT ALPHA)
		optimizeGlobalRate();
	} else if (_args_info.optimizeAlpha_flag) { // here we do alpha but NO bbl, and no global rate.
		optimizeAlphaOnly();
	}
}


void mainSemphy::compute(bool doBootstrapOneIteration) {
	// The program has layers. The first layer computs a tree.
	computeTree(doBootstrapOneIteration);

	// The second layer is parameter optimization
	optimizeParameters(); //(branch lengths, alpha, both, ...) 
	
	// The third layer computes some likelihoods on a given tree.
	computeLikelihoodAndLikelihoodPerPositionAndPosterior();
}

void mainSemphy::output() const {
	// The forth layer prints results.
	out()<<"# Finished tree reconstruction."<<endl;
	LOG(3,<<"# Finished tree reconstruction."<<endl);
	printLikelihoodAndLikelihoodPerPosition();
	printPosterior();
	printTree();
	printTreeToTreeFile();
}

// this function is NOT needed now for SEMPHY, per-ce, but is used with SEMPHY LIB.
void mainSemphy::extractPosteriorProb(VVdoubleRep & posteriorProbVV) const {
	likelihoodComputation::getPosteriorOfRates(*_etPtr,_sc,*_spP,posteriorProbVV,_weights);
}
