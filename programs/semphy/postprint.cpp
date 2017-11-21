// $Id: postprint.cpp 409 2005-06-28 13:12:24Z ninio $

#include "tree.h"
#include "sequenceContainer.h"
#include "mainSemphy.h"
#include "njGamma.h"
#include "logFile.h"
#include "treeUtil.h"
#include "likeDist.h"
#include "someUtil.h" 
#include "codon.h"
#include "recognizeFormat.h"

#include "pDistance.h"
#include "jcDistance.h"
#include "likeDist.h"


#include <iostream>
#include <cassert>
using namespace std;


int main(int argc,char* argv[]) {
  if (argc < 2) errorMsg::reportError("The program must get some parameters in the command line, use -h for help");
  semphy_args_info args_info;
  if (semphy_cmdline_parser(argc, argv, &args_info) != 0) {
    errorMsg::reportError("error reading command line",1);
  }
  
   // check that we have all we need:

  if (!args_info.tree_given) {
    errorMsg::reportError("Must enter input tree",1);
  }
 
  mainSemphy ms(args_info);
   assert(ms.getStochasticProcess().categories() ==1);
   for (MDOUBLE f=0.005;f<5;f+=0.005){
     ms.setGlobalRate(f);
     ms.out()<<"global rate "<<f<<"  "<<     ms.getStochasticProcess().getGlobalRate()<<endl;
     ms.computeLikelihoodAndLikelihoodPerPosition();
     ms.printLikelihoodAndLikelihoodPerPosition();
   }
  return 0;
}
