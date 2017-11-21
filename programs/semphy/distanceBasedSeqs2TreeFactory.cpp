// 	$Id: distanceBasedSeqs2TreeFactory.cpp 1941 2007-04-15 14:17:22Z privmane $	
#include "distanceBasedSeqs2TreeFactory.h"
#include "ssrvDistanceSeqs2Tree.h"

distanceBasedSeqs2Tree* distanceBasedSeqs2TreeFactory(const distanceBasedMethod_t distanceBasedMethod, 
													  stochasticProcess& sp, // may change sp (alpha and/or nu)
													  const bool   useJcDistance,
													  const bool   optimizeAlpha,
													  const bool   useSSRV,
													  const double epsilonLikelihoodImprovement4iterNJ,
													  const double epsilonLikelihoodImprovement4pairwiseDistance,
													  const double epsilonLikelihoodImprovement4alphaOptimiz,
													  const double epsilonLikelihoodImprovement4BBL,
													  const int maxIterationsBBL)
{

	// 1. Construct the object of the requested distance-based tree-reconstruction method (class distances2Tree)
	NJalg NJa;
	distances2Tree* d2tPtr = &NJa; // for future polymorphism

	// 2. Construct the object of the requested distance estimation method (class distanceMethod) distanceMethod *dmPtr = NULL;
	// And construct the object of the distance-based sequences-to-tree method (class distanceBasedSeqs2Tree)
	//	distanceMethod *dmPtr = NULL;
	distanceBasedSeqs2Tree *s2tPtr = NULL;
	
	
	switch (distanceBasedMethod) {
	case homogeneousRatesDTME:	
		if (useJcDistance) {
			jcDistance jcd;
			s2tPtr = new distanceBasedSeqs2Tree(jcd, *d2tPtr);  
		} else {
			likeDist ld(sp,epsilonLikelihoodImprovement4pairwiseDistance);
			s2tPtr = new distanceBasedSeqs2Tree(ld, *d2tPtr);  

		}
		break;
	case pairwiseGammaDTME:
		if (optimizeAlpha) {
			pairwiseGammaDistance pgd(sp,epsilonLikelihoodImprovement4pairwiseDistance); // distance method that optimizes alpha for the pair of sequences
			s2tPtr = new distanceBasedSeqs2Tree(pgd, *d2tPtr);
		} else {
			likeDist ld1(sp,epsilonLikelihoodImprovement4pairwiseDistance); // distance method that uses the given alpha with no optimization
			s2tPtr = new distanceBasedSeqs2Tree(ld1, *d2tPtr);
		}
		break;
	case commonAlphaDTME: {
		if (!useSSRV) {
			likeDist ld2(sp,epsilonLikelihoodImprovement4pairwiseDistance);
			s2tPtr = new commonAlphaDistanceSeqs2Tree(ld2, *d2tPtr, NULL, 
													  epsilonLikelihoodImprovement4iterNJ,
													  epsilonLikelihoodImprovement4alphaOptimiz,
													  epsilonLikelihoodImprovement4BBL,
													  maxIterationsBBL); 
		} else {
			likeDist ld3(sp,epsilonLikelihoodImprovement4pairwiseDistance);
			s2tPtr = new ssrvDistanceSeqs2Tree(ld3, *d2tPtr, NULL,
											   epsilonLikelihoodImprovement4iterNJ,
											   epsilonLikelihoodImprovement4alphaOptimiz,
											   epsilonLikelihoodImprovement4BBL,
											   maxIterationsBBL);

		}
	}
		break;
	case rate4siteDTME: {
		givenRatesMLDistance grd(sp,epsilonLikelihoodImprovement4pairwiseDistance);
		s2tPtr = new rate4siteDistanceSeqs2Tree(grd, *d2tPtr, NULL,
												epsilonLikelihoodImprovement4iterNJ,
												epsilonLikelihoodImprovement4alphaOptimiz,
												epsilonLikelihoodImprovement4BBL,
												maxIterationsBBL); }
		break;
	case posteriorDTME: {
		posteriorDistance posd(sp,epsilonLikelihoodImprovement4pairwiseDistance);
		s2tPtr = new posteriorDistanceSeqs2Tree(posd, *d2tPtr, NULL,
												epsilonLikelihoodImprovement4iterNJ,
												epsilonLikelihoodImprovement4alphaOptimiz,
												epsilonLikelihoodImprovement4BBL,
												maxIterationsBBL);}
		break;
	default: errorMsg::reportError("bad distanceBasedMethod - method not in the list of implimented methods");
	}
	return s2tPtr;
}






