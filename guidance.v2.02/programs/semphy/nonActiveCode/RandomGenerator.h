#ifndef RandomGenerator_h
#define RandomGenerator_h

// Marsaglia's subtractive R.N. generator with carry; combined with a weyl
// generator.
// Source: Computer Physics Communications 60 (1990) 345-349.
// Written by Geoffrey Zweig, 1992.
// IMPLEMENTATION FILE

#ifdef __GNUG__
#ifdef PRAGMA_TEMPLATES
#pragma interface
#endif
#endif

#include "Multinomial.h"
#include <cmath> //CHANGE
//#include "general.h"
class tMarsagliaGenerator
{
public:
  tMarsagliaGenerator(unsigned query=0);
  void Initialize(unsigned start);
  unsigned long RandomLong() {return Next();} 

  unsigned long RandomLong(unsigned long range)
  {return (Next() % range);}

  unsigned RandomInt(unsigned range) {return unsigned(Next() % range);} 

  double RandomDouble(double range);	
  
private:
  unsigned long word1[44];
  unsigned long weyl;
  int i,j,carry;
  unsigned long Next();
};


class tRandomGenerator : public tMarsagliaGenerator {
public:



  // Sample Uniform (0,1)
  double SampleUniform();

  // sample from a multinomial
  int SampleMultinomial( int n, double const* );
  int sampleMulti(tMultinomial const & m );

  vector<int> sampleGroup(int size,int groupSize);
  
  // Array containing log probabilities (the array is modified !)
  int SampleLogMultinomial( int n, double* );

  // Sample from Normal distribution

  double SampleNormal();

  double SampleGaussian(double mu, double prec);
  
  // Sample from a Gamma distribution

  double SampleGamma(double Alpha);

  double SampleGamma(double Alpha, double Beta);

  
  // Sample from a Dirichlet Distribution
  
  void   SampleDirichlet(const vector<double>& rgalpha,
			 vector<double>& rgprob);
  
  
  void   SampleDirichlet(const vector<double>& mean,
			 double precision,
			 vector<double>& rgprob );

  void   SampleDirichlet( int n,
			  double const*alpha,
			  double* prob );
  
  void   SampleDirichlet( int n,
			  double const* mean,
			  double        precision,
			  double* prob );

  void   SampleDirichlet( int n,
			  float const* mean,
			  double        precision,
			  float* prob );

protected:
  double DblGammaGreaterThanOne(double dblAlpha);
  double DblGammaLessThanOne(double dblAlpha);
  
};

extern tRandomGenerator _RandomProbGenerator;

inline
double
tRandomGenerator::SampleGaussian(double mu, double prec)
{
  double sigma = 1/sqrt(prec);
  double x= SampleNormal()*sigma + mu;

#ifdef LOG  
  cerr << "SampleGaussian(" << mu << " " << sigma << ") = " << x << "\n";
#endif
  
  return x;
}

inline
double
tRandomGenerator::SampleGamma(double Alpha, double Beta)
{
  double x= SampleGamma(Alpha)/Beta;

#ifdef LOG  
  cerr << "SampleGamma(" << Alpha << " " << Beta << ") = " << x << "\n";
#endif
  
  return x;
}
inline
void
tRandomGenerator::SampleDirichlet(const vector<double>& rgalpha,
				  vector<double>& rgprob)
{
  SampleDirichlet( rgalpha, 1.0, rgprob );
}

inline
void
tRandomGenerator::SampleDirichlet( int n,
				   double const* alpha,
				   double * prob )
{
  SampleDirichlet( n, alpha, 1.0, prob );
}


#endif


