// $Id: semphyStep.cpp 7805 2010-03-26 13:59:46Z privmane $

#include "semphyStep.h"
#include "computePijComponent.h"
#include "approxSemphyDistance.h"
#include "minSpanTree.h"
#include "rearangeTree.h"
#include "correctToCanonialForm.h"
#include "likelihoodComputation.h"
#include "searchStatus.h"
#include "talRandom.h"
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "constraints.h"
#include "treeIt.h"
#include "someUtil.h"
#include <cstdlib>

//#include "bblEM.h"
//#include "computeExactAlg.h"
//#include "computeMarginalAlg.h"
//#include "exactSemphyDistance.h"
//#include "jointTable.h"

semphyStep::semphyStep(tree& et,const tree* ctPtr,
		  const sequenceContainer& sc,
		  const stochasticProcess& sp,
		  const computePijGam& pij0,
		  const suffStatGlobalGam& cup,
		  const suffStatGlobalGam& cdown, const bool useApproxCounts,
		  const VdoubleRep& cprobAtEachPos,
		  const VVdoubleRep& posteriorRateProbAtEachPos,
		  const suffStatGlobalGam& computeMarginal,
				const Vdouble *weights,
				const MDOUBLE toll) : _et(et),_ctPtr(ctPtr),_sc(sc),_sp(sp),_pij0(pij0),_cup(cup),_cdown(cdown),
				_useApproxCounts(useApproxCounts),_cprobAtEachPos(cprobAtEachPos),
				_posteriorRateProbAtEachPos(posteriorRateProbAtEachPos),
				_computeMarginal(computeMarginal),_weights(weights),_toll(toll){
	computeSemphyStep();
}

// this function fixes the penalty table, so there is HUGE penalty
// for connecting any member of the clade to a non-member.
// We note, that the members of a "constrain clade" are all the leaves
// in this clade, but also all the internal nodes of the clade,
// including the root of the clade.
// The idea is that by doing so, all the members of the clade are sure
// to be clustered together, and all the members outside are also
// sure to be clustered together.
// Anyway, the HUGE penalty must be taken in the spanning tree
// at least once: when the clade is connected to the rest of the "world".
// Because we put a "HUGE" penalty, this penalty will only be taken once.
// This will also work in the general case of several (nested or not) clades.

// A condition that must be verify in order for the function to work:
// say we want sequence A to cluster with sequence B.
// What is the internal node connecting these two sequences?
// In order for the algorithm to work, we must demand that the tree
// according to which the penalty matrix was computed is compatible with the
// constraint tree.
// Thus, we take the internal nodes that are relevant - from the tree upon which
// the constraints were computed.

void semphyStep::addConstraintPenalty(VVdouble & penaltyTable) {
	// Here we set the constraints.
	// One of the functionality of this class, is to create a penalty
	// table (VVdouble) that corresponds to the constraints.
	// For example, if the constraint is that seq A and B are together,
	// Say that node N is their parent in the tree.
	// If S = {A,B,N}
	// In the penalty table that is created, there will be a "1"
	// between and member of S to a non member of S.
	// In all other cases, there will be 0.
	// If there are multiple constraints - it will be a superposition
	// of all the constaints. For example if
	// S1 = {A,B,C,N1,N2}
	// S2 = {A,B,N1}
	// In the penalty table between A and B it will be 0.
	// between A and C, it will be 1.
	// Between A and a node that is not in S1 or in S2 - it will be 2.
	constraints cons(*_ctPtr);
	cons.setTree(_et);

	// here we check the condition described above is fulfilled.
	if (!cons.fitsConstraints()){ // sanaty check
	  LOGDO(1,_et.output(myLog::LogFile(),tree::PHYLIP,true));
	  LOGDO(1,_ctPtr->output(myLog::LogFile(),tree::PHYLIP,true));
	  LOGDO(1,cons.outputMissingClads(cerr));
	  errorMsg::reportError(" Tree does not fit constraints in SEMPHY step");
	}

	VVdouble consPeneTable(cons.getPeneltyTable());

	// here we compute the sum of penalty
	double sumOfPenalties=0.0;
	for (int i=0;i<penaltyTable.size();++i) {
		for (int j=0;j<penaltyTable[i].size();++j) {
			sumOfPenalties+=penaltyTable[i][j];
		}
	}

	// here the motivation is to create a penalty
	// that is two order of magntide bigger than the
	// sum of penalty.
	// In order for the penalty number to be a nice integer
	// we use the following trick.
	// Consider for example the case where the sum was 7
	// Log10(7) = 0.845;
	// Then, log10(7) + 2 = 2.845
	// ceil just rounds up to the next integer
	// ceil(2.845) = 3
	// Then pow(10.0,3) = 1000.
	// which is a nice number, about 2+ order of magnitude
	// bigger than 7.
	double LargePenelty =pow(10.0,ceil(log10(sumOfPenalties)+2.0));


	for (int i2=0;i2<penaltyTable.size();++i2)
	  for (int j=0;j<penaltyTable[i2].size();++j)
	    penaltyTable[i2][j]-=consPeneTable[i2][j]*LargePenelty;
}

void semphyStep::computeSemphyStep(){
	// for the starting tree decide which nodes are "real" and which are not.
	vector<char> isRealTaxa(_et.getNodesNum(),0);
	treeIterTopDownConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf()) {isRealTaxa[mynode->id()] = 1;}
	}

//	fillCpijUpDownExactMarginalProbAtEachPos();
	semphyDistance* semDis1 = computeQmatrix(_weights);// NULL = WEIGHTS...
	VVdouble penetlyTable = *(semDis1->getLikeDistanceTablePtr()); // copy
	MDOUBLE score = 0.0;
	if (_ctPtr!=NULL) {		
	  addConstraintPenalty(penetlyTable);
	}
	rearrangeTree::pairSet inSet = minSpanTree::span_tree(penetlyTable,&score);

	maybePrintPenaltyTableAndStartingTree(penetlyTable); //(log output lvl 50)
	MDOUBLE qSTART = computeQstart(semDis1);
	maybePrintDistanceBetweenNodesTable(semDis1);//(log, 50)
	maybePrintSpanTreeListAndQAfterSpanTree(inSet,score);//(log, 50)

	rearrangeTree rearrangeTree1(&inSet,semDis1);
	rearrangeTree1.reconstructTree(_et);
	maybePrintTheTreeAfterRearrangeTree();//(log, 50)
	correctToCanonialForm ctcf(&_et,*semDis1->getDistanceTablePtr(),isRealTaxa);
	ctcf.correctTree();
	maybePrintTheTreeAfterCorrectToCanonialForm();//(log, 50)
	maybePrintQspanMinusQinit(score,qSTART);//(log, 50)
	if (semDis1!=NULL) delete semDis1;
}
/*
void semphyStep::fillCpijUpDownExactMarginalProbAtEachPos() {
	computeUpAlg::fillComputeUp(&_et,_pi,_computeUp1);
	_computeProbOfEachPos1->fillProbOfEachPosition(&_et,_pi,_computeUp1);
	computeDownAlg::fillDown(&_et,_pi,_computeUp1,_computeDown1);
	computeExactAlg::fillExact(&_et,_pi,_computeUp1,_computeDown1,_computeExact1);
	computeMarginalAlg::fillMarginal(&_et,_pi,_computeExact1,_computeProbOfEachPos1,_computeMarginal1);
}*/

void semphyStep::maybePrintPenaltyTableAndStartingTree(const VVdouble& penetlyTable) const{
	LOG(50,<< "start of semphy step"<<endl);
	LOGDO(50,_et.output(myLog::LogFile(),tree::ANCESTOR));
	LOGDO(50,_et.output(myLog::LogFile(),tree::PHYLIP));
	LOG(50, <<"START L = "<<likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,NULL)<<endl);
	vector<tree::nodeP> allNodes;
	_et.getAllNodes(allNodes,_et.getRoot());
	for (int n1=0; n1 < allNodes.size(); ++n1) 
		LOG(50,<<allNodes[n1]->name()<<"\tid= "<<allNodes[n1]->id()<<endl);
	LOG(50,<<"the penetlyTable: "<<endl);
	for (int i=0; i < penetlyTable.size(); ++i) {
		for (int k=0; k < penetlyTable[i].size(); ++k) {
			LOG(50,<< penetlyTable[i][k]<<"\t");
		}
		LOG(50,<<endl);
	}
}

MDOUBLE semphyStep::computeQstart(const semphyDistance* inSemDist) const {
	MDOUBLE qSTART=0.0;
	vector<tree::nodeP> allNodes;
	_et.getAllNodes(allNodes,_et.getRoot());
	for (int n1=0; n1 < allNodes.size(); ++n1) {
		if (allNodes[n1]->father()!=NULL) qSTART += inSemDist->getLikeDistance(allNodes[n1]->id(),allNodes[n1]->father()->id());
	}
	LOG(50,<<" q start = "<<qSTART<<endl);
	return qSTART;
}
	
void semphyStep::maybePrintDistanceBetweenNodesTable(const semphyDistance* inSemDist) const  {
	if (50<=myLog::LogLevel()) {
		LOG(50,<< "printing the branch length table"<<endl);
		VVdouble tTable = *inSemDist->getDistanceTablePtr(); // copy
		for (int i=0; i < tTable.size(); ++i) {
			for (int k=0; k < tTable[i].size(); ++k) {
				LOG(50,<< tTable[i][k]<<"\t");
			}
			LOG(50,<<endl);
		}
	}
}

void semphyStep::maybePrintSpanTreeListAndQAfterSpanTree(
		const rearrangeTree::pairSet& inSet,const MDOUBLE score) const{
	if (50<=myLog::LogLevel()) {
		LOG(50,<<"printing span tree connection (id) set"<<endl);
		for (rearrangeTree::pairSet::const_iterator z= inSet.begin(); z!= inSet.end(); ++z) {
			LOG(50,<<(z->first)<<" is connected to "<< z->second<<endl);
		}
		LOG(50,<<endl);
		LOG(50,<<"the score of the span tree: "<<score<<endl);
	}
}	

void semphyStep::maybePrintTheTreeAfterRearrangeTree() const{
	LOG(50,<<"after rearrange tree"<<endl);
	LOGDO(50,_et.output(myLog::LogFile(),tree::ANCESTOR));
	LOGDO(50,_et.output(myLog::LogFile(),tree::PHYLIP));
}

void semphyStep::maybePrintTheTreeAfterCorrectToCanonialForm() const{
	LOG(50,<<"after correction to canonial form"<<endl);
	LOGDO(50,_et.output(myLog::LogFile(),tree::ANCESTOR));
	LOGDO(50,_et.output(myLog::LogFile(),tree::PHYLIP));
}

void semphyStep::maybePrintQspanMinusQinit(const MDOUBLE qpan, const MDOUBLE qstart) const{
	LOG(50,<<" q span - q start = "<<qpan-qstart<<endl);
}

semphyDistance* semphyStep::computeQmatrix(const Vdouble* weight) {
//	_searchStat1->createRandomWeightsGamma(_t1->seqLen());
	semphyDistance* semDis1 = NULL;
	if (_useApproxCounts) {
		//cerr<<" using APPROX! "<<endl;
		semDis1 =  new approxSemphyDistance(_et,
						_sc,
						_sp,
						_pij0,
						_cup,
						_cdown,
						_cprobAtEachPos,
						_posteriorRateProbAtEachPos,
						_computeMarginal,
						_weights,
						_toll);
	}
	else {
		cerr<<" using EXACT!. To implement again... "<<endl;
		exit(4);
	  //semDis1 =  new exactSemphyDistance(_et,
		//			     *_pi,
		//			     *_computeUp1,
		//			     *_computeDown1,
		//			     *_computeProbOfEachPos1,
		//			     *_computeMarginal1,
		//			     weight,
		//				 _toll);
	}
	semDis1->computeDistances();
	return semDis1;
}

