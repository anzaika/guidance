// $Id: njGammaMain.cpp 674 2006-04-09 17:06:42Z ninio $

#include "tree.h"
#include "sequenceContainer.h"
#include "mainSemphy.h"
#include "njGamma.h"
#include "logFile.h"
#include "treeUtil.h"
#include "likeDist.h"
#include "someUtil.h"

#include <iostream>
#include <cassert>
using namespace std;


int main(int argc,char* argv[]) {
  bool doPairwiseGammaMLDis = false; // do posterior.
// bool doPairwiseGammaMLDis = true;
  bool withoutBBL = false;
  if (argc < 2) errorMsg::reportError("The program must get some parameters in the command line, use -h for help");
  semphy_args_info args_info;
  if (semphy_cmdline_parser(argc, argv, &args_info) != 0) {
    errorMsg::reportError("error reading command line",1);
  }

  // we need a starting tree to get starting posterior probs.
  doPairwiseGammaMLDis=args_info.ADVNoPost_flag;
  args_info.optimizeAlpha_flag = 1; // true
  args_info.alpha_given = 1;
  mainSemphy ms(args_info);
  assert(ms.getStochasticProcess().categories() >1);
  ms.computeNJtree();
  if (withoutBBL == false) {
    ms.optimizeBranchLengths(); // this will do both bbl and alpha.
  } else {
    ms.optimizeAlphaOnly();
  }
  //  ms._args_info.optimizeGamma_flag = 0; // stop doing alpha optimization.
	
  VVdouble posteriorProbVV; // pos * rate

  // now we need a distance table:
  MDOUBLE newL = VERYSMALL;
  MDOUBLE oldL = VERYSMALL;
  tree oldT = ms.getTree();
  tree newT = ms.getTree();
  const int maxNumIter = 50;
  bool converged = false;
  int iterNum=0;
  distanceMethod * disMethod;
  while ( !converged && iterNum<maxNumIter) {
    stochasticProcess tmpSP = ms.getStochasticProcess();
    if (doPairwiseGammaMLDis==false) {
      LOG(5,<<"Gamma-NJ-POSTERIOR"<<endl);
      ms.extractPosteriorProb(posteriorProbVV);
      LOG(10,<<posteriorProbVV<<endl);
      disMethod = new gammaMLDistances(tmpSP,posteriorProbVV);
    } else {
      LOG(5,<<"Gamma-NJ-PAIRWISE"<<endl);
      disMethod = new likeDist(tmpSP,0.001);
    }
    VVdouble disTab;
    vector<string> vNames;
    ms.getDistanceTableAndNames(disTab,vNames,disMethod);
    
    LOG(10,<<"Distance table"<<endl<<disTab<<endl<<endl);
    for (int vNames_I=0;vNames_I<vNames.size();++vNames_I){
      LOG(10,<<vNames[vNames_I]<<endl);
    }
    LOG(10,<<endl<<endl);
      
    ms.computeNJtreeFromDisTableAndNames(disTab,vNames); 
    newT = ms.getTree();
		
    if (sameTreeTolopogy(newT,oldT)) {
      converged = true;
    } else {
      ms.computeLikelihoodAndLikelihoodPerPosition();
      newL = ms.getLikelihood();
      LOG(5,<<"Iteration = "<<iterNum<<" LL = "<<newL<<endl);
      LOG(5,<<"Found a new tree topology!"<<endl);
      LOGDO(5,newT.output(myLog::LogFile()));
      if (withoutBBL==false) {
	LOG(10,<<"optimize gamma? "<<ms._args_info.optimizeGamma_flag<<endl);
	ms.optimizeBranchLengths(); // this will do bbl 
	LOG(5,<<"Iteration = "<<iterNum<<" LL after bbl = "<<newL<<endl);
	newL = ms.getLikelihood();
	if (newL > oldL) {
	  oldL = newL;
	  oldT = newT;
	} else {
	  LOG(5,<<"Stoping GAMMA-NJ Because LL decreased, though topology has changed!");
	  converged = true;
	  ms.setTree(oldT);
	}
      } else {
	oldT = newT;
      }
    }
    iterNum++;
    delete disMethod;
  }
  if (iterNum==maxNumIter) {
    LOG(3,<<"Stoped iterating not because of convergence, but because max num of iterations");
  }

  ms.optimizeBranchLengths();
  newL = ms.getLikelihood();
  LOG(5,<<"Iteration = "<<iterNum<<" FINAL LOG LIKELIHOOD = "<<newL<<endl);

  // print the likelihood and the resulting tree and alpha.
  ms.output();
  return 0;
}
