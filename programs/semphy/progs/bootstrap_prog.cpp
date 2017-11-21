#include "bootstrap.h"
#include "treeUtil.h"

#include "constantSemphyDistance.h"
#include "rearangeTree.h"
#include "logFile.h"
#include "bootstrap_prog_cmdline.h"


void setLog(const char* logfilename, const int loglvl)
{
  if (!strcmp(logfilename,"-")||!strcmp(logfilename,""))
	  {
	    myLog::setLogOstream(&cout);
	  }
	else
	  {
	    ofstream* outLF = new ofstream;
	    outLF->open(logfilename);
	    if (!outLF->is_open()) {
	      errorMsg::reportError("unable to open file for reading");
	    }
	    myLog::setLogOstream(outLF);
	  }
	myLog::setLogLvl(loglvl);
	LOG(3,<<"START OF MUL LOG FILE\n\n");
}



int main(int argc, char *argv[])
{
  // STEP 2: WE BUILD THE STRUCT OF SEMPHY INFO
  bootstrap_prog_args_info args_info;
  if (bootstrap_prog_cmdline_parser (argc, argv, &args_info) != 0) {
    errorMsg::reportError("error in command line parsing",1);
  }
  setLog(args_info.Logfile_arg, args_info.verbose_arg);


  vector<tree> tv(getStartingTreeVecFromFile(args_info.treesList_arg));
  bootstrap b1(tv);

  if (args_info.reftree_given) {
    tree reftree(args_info.reftree_arg);
    map<int, MDOUBLE> m1(b1.getWeightsForTree(reftree)) ;
    b1.printTreeWithBPvalues(cout,reftree,m1,!args_info.noBranchLen_flag);
    cout <<endl;
  } else {			// constract tree from list
    tree t(b1.consensusTree(args_info.ConsensusLevel_arg));

	map<int, MDOUBLE> m2(b1.getWeightsForTree(t));
						 
    b1.printTreeWithBPvalues(cout,t,m2,!args_info.noBranchLen_flag);
    cout <<endl;

//     tree t(tv[0]);
  
//     semphyDistance* semDis1 = NULL;
//     semDis1 = new constantSemphyDistance();

//     set< pair<int,int> > outset;
//     vector<float> support=b1.thresholdTree(t,outset, args_info.ConsensusLevel_arg);
//     rearrangeTree  rt(&outset,semDis1);
//     rt.reconstructTree(t);
//    b1.printTreeNH(cout, t, support);
  }
  return (0);
}
