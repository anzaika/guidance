// $Id: posterior.cpp 674 2006-04-09 17:06:42Z ninio $

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
  bool withoutBBL = false;
  if (argc < 2) errorMsg::reportError("The program must get some parameters in the command line, use -h for help");
  semphy_args_info args_info;
  if (semphy_cmdline_parser(argc, argv, &args_info) != 0) {
    errorMsg::reportError("error reading command line",1);
  }


  // check that we have all we need:

  if (!args_info.tree_given) {
    errorMsg::reportError("Must enter input tree",1);
  }
  if (!args_info.alpha_given) {
    errorMsg::reportError("Must set gamma",1);
  }



  mainSemphy ms(args_info);
  // see that we actually got gamma
  assert(ms.getStochasticProcess().categories() >1);
  for (int i=0;i<ms.getStochasticProcess().categories();++i)
    LOG(3,<<"category["<<i<<"]="<<ms.getStochasticProcess().rates(i)<<endl);
  if (withoutBBL == false) {
    ms.optimizeBranchLengths(); // this will do both bbl and alpha.
  } 
	
  VVdouble posteriorProbVV; // pos * rate
  ms.extractPosteriorProb(posteriorProbVV);
  ms.out()<<posteriorProbVV<<endl;
  return 0;
}
