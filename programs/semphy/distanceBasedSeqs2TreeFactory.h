// 	$Id: distanceBasedSeqs2TreeFactory.h 1941 2007-04-15 14:17:22Z privmane $	
#include "distanceBasedSeqs2Tree.h"

#include "pairwiseGammaDistance.h"

#include "nj.h"

//using namespace std;

#ifndef __DISTANCEBASEDSEQS2TREEFACTORY_H
#define __DISTANCEBASEDSEQS2TREEFACTORY_H

typedef enum {homogeneousRatesDTME,
			  pairwiseGammaDTME,
			  commonAlphaDTME,
			  rate4siteDTME,
			  posteriorDTME
} distanceBasedMethod_t;



distanceBasedSeqs2Tree* distanceBasedSeqs2TreeFactory(const distanceBasedMethod_t distanceBasedMethod, 
													  stochasticProcess& sp, // may change sp (alpha)
													  const bool   useJcDistance,
													  const bool   optimizeAlpha,
													  const bool   useSSRV,
													  const double epsilonLikelihoodImprovement4iterNJ = 0.01,
													  const double epsilonLikelihoodImprovement4pairwiseDistance = 0.0001,
													  const double epsilonLikelihoodImprovement4alphaOptimiz = 0.01,
													  const double epsilonLikelihoodImprovement4BBL = 0.01,
													  const int maxIterationsBBL = 10);


#endif





