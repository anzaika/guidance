// $Id: semphy.cpp 1282 2006-12-05 11:03:10Z privmane $

#include "bootstrapProvider.h"
#include "mainSemphy.h"
#include "semphy_cmdline.h"
#include "errorMsg.h"

int main(int argc, char* argv[]) {
	if (argc < 2) errorMsg::reportError("The program must get some parameters in the command line, use -h for help");
	gengetopt_args_info args_info;
	if (cmdline_parser(argc, argv, &args_info) != 0) {
	  errorMsg::reportError("error reading command line",1);
	}

	if (args_info.consurf_flag){
	  tree t_tmp;
	  {			// for ms1;
	    args_info.optimizeAlpha_flag = 0; // true
	    args_info.alpha_given = 0;
	    mainSemphy ms_tmp(args_info);
	    ms_tmp.computeNJtree();
	    t_tmp=ms_tmp.getTree();
	  }
	  args_info.optimizeAlpha_flag = 1; // true
	  args_info.alpha_given = 1;
  	  mainSemphy ms2(args_info);
	  myLog::printArgv(1, argc, argv);
  	  ms2.setTree(t_tmp);
  	  ms2.optimizeAlphaOnly();
  	  args_info.alpha_arg =   static_cast<gammaDistribution*>(ms2.getStochasticProcess().distr())->getAlpha();
	}
	mainSemphy ms(args_info);
	myLog::printArgv(1, argc, argv);
	if (!args_info.BPonUserTree_given) {
		// compute the ML tree for example.
		ms.compute();
		ms.output();
	}

	// computing the BP values
	if (args_info.BPrepeats_given) {
		bootstrapProvider bp(args_info);
		bp.computeBP(ms);
		bp.output(ms.out());
	}
	return 0;
}
