// $Id: distances.cpp 674 2006-04-09 17:06:42Z ninio $

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
  
  if (args_info.inputs_num<1) 
    errorMsg::reportError("Need a file with an extra pair of sequances",1);



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
    LOG(3,<<"category["<<i<<"]="<<ms.getStochasticProcess().rates(i)<<"\t"<<ms.getStochasticProcess().ratesProb(i)<<endl);
  if (args_info.likelihood_flag)
    ms.computeLikelihoodAndLikelihoodPerPosition();
  else
    ms.optimizeBranchLengths(); // this will do both bbl and alpha.

  ms.printLikelihoodAndLikelihoodPerPosition();

  ms.printTree();

  VVdouble posteriorProbVV; // pos * rate
  ms.extractPosteriorProb(posteriorProbVV);
  LOGDO(3,  ms.out()<<posteriorProbVV<<endl);

  // read the extra sequances
  sequenceContainer sc;
  alphabet* alphP=NULL;
  switch (args_info.alphabet_arg) 
    { // allwayes defined, with default
    case 4:	
      alphP = new nucleotide; 
      break;
    case 20: 
      alphP = new amino; 
      break;
    case 64: 
      alphP = new codon; 
	    break;
    default: errorMsg::reportError("error getting the alphabet");
    }

  if (args_info.categories_given)
    LOG(3,<<"using "<<args_info.categories_arg<<" bins for gamma"<<endl); 
  if (args_info.inputs_num <1)
    errorMsg::reportError("please provide a file with a list of sequnce pair files");
  ifstream seqlist;
  seqlist.open(args_info.inputs[0]);
  if (! seqlist.good())
    errorMsg::reportError("can not open sequence pair list file");
  char tmp[1000];		// this is not ideal!

  while (!seqlist.eof()){
    seqlist >> tmp;
    if (seqlist.eof()) break;
    string seqid(tmp);
    ms.out() << seqid <<"\t";

    seqlist >> tmp;
    string sequenceFileName(tmp);
    ms.out() << sequenceFileName<<"\t";

    ifstream ins;
    ins.open(sequenceFileName.c_str());
    if (! ins.good())
      errorMsg::reportError("can not open sequence pair file");
    
    sc = recognizeFormat::read(ins,alphP);
    
    ins.close();
      
    
    
    jcDistance jcD(args_info.alphabet_arg);
    ms.out() <<""<<jcD.giveDistance(sc[0], sc[1], NULL)<<"\t";
    
    pDistance pD;
    ms.out() <<""<<pD.giveDistance(sc[0], sc[1], NULL)<<"\t";
    
    pupAll probMod(datMatrixHolder::dayhoff);
    chebyshevAccelerator pijAcc(&probMod);


    uniDistribution dist1;
    stochasticProcess sp1(&dist1, &pijAcc);
    likeDist lD1(sp1);
    ms.out() <<""<<lD1.giveDistance(sc[0], sc[1], NULL)<<"\t";

    gammaDistribution dist2(args_info.gamma_arg,args_info.categories_arg);
    stochasticProcess sp2 = ms.getStochasticProcess();
    likeDist lD2(sp2);
    ms.out() <<""<<lD2.giveDistance(sc[0], sc[1], NULL)<<"\t";

    gammaMLDistances lD3(sp2,posteriorProbVV);
    ms.out() <<""<<lD3.giveDistance(sc[0], sc[1], NULL)<<endl;

  }

  delete alphP;
  return 0;
}
