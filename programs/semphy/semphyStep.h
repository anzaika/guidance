// $Id: semphyStep.h 702 2006-05-27 17:16:58Z ninio $

#ifndef ___SEMPHY_STEP
#define ___SEMPHY_STEP

#include "suffStatComponent.h"
#include "semphyDistance.h"
#include "rearangeTree.h"
#include "minSpanTree.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"
#include "stochasticProcess.h"

class semphyStep {
public:
	explicit semphyStep(
		  tree& et,
		  const tree* ctPtr,
		  const sequenceContainer& sc,
		  const stochasticProcess& sp,
		  const computePijGam& pij0,
		  const suffStatGlobalGam& cup,
		  const suffStatGlobalGam& cdown,
		  const bool useApproxCounts,
		  const VdoubleRep& cprobAtEachPos,
		  const VVdoubleRep & posteriorRateProbAtEachPos,
		  const suffStatGlobalGam& computeMarginal,
				const Vdouble *weights,
				const MDOUBLE toll);

  void doNothing(void) const {};
				      

private:

	void computeSemphyStep();
//	void fillCpijUpDownExactMarginalProbAtEachPos();
	semphyDistance* computeQmatrix(const Vdouble* weight);
	MDOUBLE computeQstart(const semphyDistance* inSemDist) const;
	void addConstraintPenalty(VVdouble & penaltyTable);
	tree& _et;
	const tree* _ctPtr;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	const computePijGam& _pij0;
	const suffStatGlobalGam& _cup;
	const suffStatGlobalGam& _cdown;
	const bool _useApproxCounts;
	const VdoubleRep& _cprobAtEachPos;
	const VVdoubleRep & _posteriorRateProbAtEachPos;
	const suffStatGlobalGam& _computeMarginal;
	const Vdouble *_weights;
	const MDOUBLE _toll;

	// ouput to log file functions...
	void maybePrintPenaltyTableAndStartingTree(const VVdouble& penetlyTable) const;
	void maybePrintDistanceBetweenNodesTable(const semphyDistance* inSemDist) const;
	void maybePrintSpanTreeListAndQAfterSpanTree(const rearrangeTree::pairSet& inSet,const MDOUBLE score) const;
	void maybePrintTheTreeAfterRearrangeTree() const;
	void maybePrintQspanMinusQinit(const MDOUBLE qpan, const MDOUBLE qstart) const;
	void maybePrintTheTreeAfterCorrectToCanonialForm() const;

};

#endif

