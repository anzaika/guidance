// $Id: semphySearchBestTree.cpp 6002 2009-03-20 19:39:03Z privmane $

#include "semphySearchBestTree.h"
#include "bestAlpha.h"
#include "uniDistribution.h"
#include "aaJC.h"
#include "nucJC.h"
#include "pijAccelerator.h"
#include "trivialAccelerator.h"
#include "chebyshevAccelerator.h"
#include "stochasticProcess.h"
#include "nucleotide.h"
#include "amino.h"
#include "sequenceContainer.h"
#include "phylipFormat.h"
#include "maseFormat.h"
#include "fastaFormat.h"
#include "tree.h"
#include "likelihoodComputation.h"
#include "computePijComponent.h"
#include "seqContainerTreeMap.h"
#include "gammaDistribution.h"
#include "bblEM.h"
#include "NNiProp.h"
#include "NNiSep.h"
#include "jcDistance.h"
#include "distanceTable.h"
#include "nj.h"
#include "semphyStep.h"
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeMarginalAlg.h"
#include "logFile.h"
#include "getRandomWeights.h"
#include "treeUtil.h"
#include "bestTamura92param.h"
#include "bestGtrModelParams.h"

using namespace std;
using namespace likelihoodComputation;


//HERE IS THE MAIN CODE FOR THE SEMPHY SEARCH
//THINGS TO TAKE INTO ACCOUNT:
//(1) Several starting point
//(2) Tmp handling for semphy-ratchet and semphy-anneal.

semphySearchBestTree::semphySearchBestTree(sequenceContainer& sc,
								  tree& startTree,
								  const tree* constraintTree,
								  stochasticProcess& sp,
								  ostream& out,
								  const int numOfRandomStart,
								  const bool optimizeAlpha,
								  const double epsilonLikelihoodImprovement4alphaOptimiz,
								  const double epsilonLikelihoodImprovement4BBL,
								  const int maxIterationsBBL,
								  const Vdouble * weights) {
	if (numOfRandomStart == 1) {
	  semphyBasicSearchBestTree(sc,startTree,constraintTree,sp,out, optimizeAlpha,
								epsilonLikelihoodImprovement4alphaOptimiz,epsilonLikelihoodImprovement4BBL,maxIterationsBBL,
								weights);
	} else {
	  semphyBasicSearchBestTreeManyRandomStarts(sc,startTree,constraintTree,sp,out,numOfRandomStart, optimizeAlpha,
												epsilonLikelihoodImprovement4alphaOptimiz,epsilonLikelihoodImprovement4BBL,maxIterationsBBL,
												weights);
	}
}



MDOUBLE semphySearchBestTree::semphyBasicSearchBestTreeManyRandomStarts(
								  sequenceContainer& sc,
								  tree& et,
								  const tree* ctPtr,
								  stochasticProcess& sp,
								  ostream& out,
								  const int nRanStart,
								  const bool optimizeAlpha,
								  const double epsilonLikelihoodImprovement4alphaOptimiz,
								  const double epsilonLikelihoodImprovement4BBL,
								  const int maxIterationsBBL,
								  const Vdouble * weights) {
  MDOUBLE bestL = semphyBasicSearchBestTree(sc,et,ctPtr,sp,optimizeAlpha,
								   epsilonLikelihoodImprovement4alphaOptimiz,epsilonLikelihoodImprovement4BBL,maxIterationsBBL,
								   weights);
	tree bestT = et;

	const MDOUBLE tmpForStartingTreeSearch = 0.1; // PUT IN THE OPTION FILE IF PROVE TO BE USEFULL...
	for (int it=0; it < nRanStart; ++it) {
	  LOG(3,<<"\n\n\n================================================================================================================\n\n\n"<<endl);
	  LOG(3,<<"best score so far: "<<bestL<<endl);
	  LOG(3,<<"starting random start number: "<<it+1<<endl);
		Vdouble startingTreeWeights(sc.seqLen(),0.0);
 		getRandomWeights::randomWeightsGamma(startingTreeWeights,tmpForStartingTreeSearch);
	
	
		jcDistance pd1; // assumes all alphabets are the same size
		VVdouble disTab;
		vector<string> vNames;
		giveDistanceTable(&pd1,
						   sc,
						   disTab,
						   vNames,
						   &startingTreeWeights);
		NJalg nj1;
		et = nj1.computeTree(disTab,vNames);
		LOG(3,<<"starting random tree for this semphy iteration "<<endl);
		LOGDO(3,et.output(myLog::LogFile()));
		MDOUBLE newL = semphyBasicSearchBestTree(sc,et,ctPtr,sp,optimizeAlpha,
								   epsilonLikelihoodImprovement4alphaOptimiz,epsilonLikelihoodImprovement4BBL,maxIterationsBBL,
								   weights);
		LOG(3,<<" with this random tree the log likelihood is: "<<newL<<endl);

		if (newL > bestL) {
		  LOG(3,<<" new best tree: score = "<<newL<<endl);

			bestT = et;
			bestL = newL;
			LOGDO(3,bestT.output(myLog::LogFile()));
		}
	}
	out<<"Semphy best tree: score = "<<bestL<<endl;
	bestT.output(out);
	return bestL;
}

MDOUBLE semphySearchBestTree::semphyBasicSearchBestTree(
								  sequenceContainer& sc,
								  tree& et,
								  const tree* ctPtr,
								  stochasticProcess& sp,
								  const bool optimizeAlpha,
								  const double epsilonLikelihoodImprovement4alphaOptimiz,
								  const double epsilonLikelihoodImprovement4BBL,
								  const int maxIterationsBBL,
								  const Vdouble * weights) {
  
  return semphyBasicSearchBestTree(sc,et,ctPtr,sp,cout,optimizeAlpha,
								   epsilonLikelihoodImprovement4alphaOptimiz,epsilonLikelihoodImprovement4BBL,maxIterationsBBL,
								   weights);
}

// Local helper function to run optimization of alpha (and other model parameters) and BBL
MDOUBLE optimizeParamAndBBL(
	sequenceContainer& sc,
	tree& et,
	stochasticProcess& sp,
	const double epsilonLikelihoodImprovement4alphaOptimiz,
	const double epsilonLikelihoodImprovement4BBL,
	const int maxIterationsBBL,
	const Vdouble * weights) 
{
	if (dynamic_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel())) {
		// optimizing params of the tamura92 model
		bestTamura92ParamAlphaAndBBL optimizer(et, sc, sp, weights, 5, 0.05,
											   epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/,
											   epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/,
											   epsilonLikelihoodImprovement4alphaOptimiz/*0.01*/,
											   epsilonLikelihoodImprovement4BBL/*0.01*/,
											   5.0, maxIterationsBBL, 1.5 );
		return optimizer.getBestL();

	} else if (dynamic_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel())) {
		// Optimizing params of the gtr model
		bestGtrModel optimizer(et, sc, sp, weights, 5,
							   epsilonLikelihoodImprovement4alphaOptimiz,
							   epsilonLikelihoodImprovement4alphaOptimiz,
							   true, true);
		return optimizer.getBestL();

	} else {
		bestAlphaAndBBL tmpbestAlpha(et, sc, sp, weights, 1.5, 5.0,
									 epsilonLikelihoodImprovement4alphaOptimiz,
									 epsilonLikelihoodImprovement4BBL,
									 maxIterationsBBL);
		return tmpbestAlpha.getBestL();
	}
}

MDOUBLE semphySearchBestTree::semphyBasicSearchBestTree(
								  sequenceContainer& sc,
								  tree& et,
								  const tree* ctPtr,
								  stochasticProcess& sp,
								  ostream& out,
								  const bool optimizeAlpha,
								  const double epsilonLikelihoodImprovement4alphaOptimiz,
								  const double epsilonLikelihoodImprovement4BBL,
								  const int maxIterationsBBL,
								  const Vdouble * weights) {
  // this is the most basic semphy search 
  //  - no anneal
  //	- no multiple random starts
  // start by optimizing the br-len of the input tree.
	MDOUBLE initial_likelihood;
  if (optimizeAlpha){
	initial_likelihood = optimizeParamAndBBL(sc, et, sp, 
											 epsilonLikelihoodImprovement4alphaOptimiz,
											 epsilonLikelihoodImprovement4BBL,
											 maxIterationsBBL, weights);
  }
  else {
    bblEM bblEM1(et,sc,sp,weights);
    //    bblEM1.getTreeLikelihood();
  }
	for (int z=0; z < 1000; ++z) {
		tree bestTree = et;
		LOG(3,<<"semphy iteration: "<<z<<endl);

		computePijGam pij0;
		pij0.fillPij(et,sp);

		suffStatGlobalGam cup;
		suffStatGlobalGam cdown;
		suffStatGlobalGam cmarg;

		cup.allocatePlace(sc.seqLen(),sp.categories(),et.getNodesNum(),sc.alphabetSize());
		cdown.allocatePlace(sc.seqLen(),sp.categories(),et.getNodesNum(),sc.alphabetSize());
		cmarg.allocatePlace(sc.seqLen(),sp.categories(),et.getNodesNum(),sc.alphabetSize());

		computeUpAlg cupAlg;
		computeDownAlg cdownAlg;
		computeMarginalAlg cmargAlg;
		int pos = 0;
		int categor = 0;

		for (pos = 0; pos < sc.seqLen(); ++pos) {
			for (categor = 0; categor < sp.categories(); ++categor) {
				cupAlg.fillComputeUp(et,sc,pos,pij0[categor],cup[pos][categor]);
				cdownAlg.fillComputeDown(et,sc,pos,pij0[categor],cdown[pos][categor],cup[pos][categor]);
			}
		}

		LOG(3,<<"starting semhy step by bbl...");
		

		VdoubleRep posLike;
		
		doubleRep tLike = likelihoodComputation::getTreeLikelihoodFromUp2(et,
							sc,
							sp,
							cup,
							posLike, // fill this vector with each position likelihood but without the weights.
							weights);//const Vdouble * weights=0);

		
		VVdoubleRep ratePosteriorProb;
		doubleRep tLike2 = likelihoodComputation::getPosteriorOfRates(et,
							sc,
							sp,
							cup,
							ratePosteriorProb, // fill this vector with each position likelihood but without the weights.
							weights);//const Vdouble * weights=0);

		// commented out by matan, 4-2-03
// 		{ // very checking
// 			for (int k=0; k < ratePosteriorProb[4].size(); ++k ) {
// 				cerr<<"ratePosteriorProb["<<k<<"]= "<<ratePosteriorProb[4][k]<<endl;
// 			}
//		}

		
		if (tLike != tLike2) { 
			cerr<<" tLike = "<<tLike<<endl;
			cerr<<" tLike2 = "<<tLike2<<endl;
			errorMsg::reportError(" error in computational - some of posterior != prob... ");
		}

		LOG(3,<<"starting likelihood = "<<tLike<<endl);
		doubleRep oldL = tLike;

		for (pos = 0; pos < sc.seqLen(); ++pos) {
			doubleRep posPr = posLike[pos];
			for (categor = 0; categor < sp.categories(); ++categor) {
				cmargAlg.fillComputeMarginal(et,sc,sp,pos,pij0[categor],cmarg[pos][categor],cup[pos][categor],cdown[pos][categor],posPr);
			}
		}



		semphyStep semphyStep1(
			  et,//tree& et,
			  ctPtr,// const tree* ctPtr,
			  sc,//const sequenceContainer& sc,
			  sp,//const stochasticProcess& sp,
			  pij0,//const computePijGam& pij0,
			  cup,//const suffStatGlobalGam& cup,
			  cdown,//const suffStatGlobalGam& cdown,
			  true,//const bool useApproxCounts,
			  posLike,//const vector<MDOUBLE>& cprobAtEachPos,
			  ratePosteriorProb, //CODE_RED
			  cmarg,//const suffStatGlobalGam& computeMarginal,
				weights,//		const Vdouble *weights,
				0.001);//	const MDOUBLE toll);

		semphyStep1.doNothing(); // just so the compiler will not complain

	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx finding the likelihood of the final trees, and output		


		MDOUBLE resx=0.0;
		if (optimizeAlpha){
			resx= optimizeParamAndBBL(sc, et, sp, 
									  epsilonLikelihoodImprovement4alphaOptimiz,
									  epsilonLikelihoodImprovement4BBL,
									  maxIterationsBBL, weights);
		}
		else{
		  bblEM bblEM1(et,sc,sp,weights);
		  resx=bblEM1.getTreeLikelihood();
		}
		 
		LOG(3,<<"likelihood after semphy step is: "<<resx<<endl);
		
		LOG(3,<<"sep best trees under semphy (gam)"<<endl);
		LOGDO(3,et.output(myLog::LogFile()));
		
		bool sameT = sameTreeTolopogy(bestTree,et);

		if (sameT) break;
		if (resx < oldL + 0.05) break;
	}

	LOG(3,<<" ============== final tree and likelihoods ================== "<<endl);
	out<<" ============== final tree and likelihoods ================== "<<endl;
	MDOUBLE resF;
	if (optimizeAlpha){
		resF= optimizeParamAndBBL(sc, et, sp, 
								  epsilonLikelihoodImprovement4alphaOptimiz,
								  epsilonLikelihoodImprovement4BBL,
								  maxIterationsBBL, weights);
	}
	else
	  {
	    bblEM bblEM1(et,sc,sp,weights);
	    resF=bblEM1.getTreeLikelihood();
	  }

		LOG(3,<<"likelihood after semphy: "<<resF<<endl);
		out<<"likelihood after semphy: "<<resF<<endl;
		
		LOG(3,<<"Best trees under semphy (gam)"<<endl);
		out<<"Best trees under semphy (gam)"<<endl;
		LOGDO(3,et.output(myLog::LogFile()));
		et.output(out);
		return resF;
}
