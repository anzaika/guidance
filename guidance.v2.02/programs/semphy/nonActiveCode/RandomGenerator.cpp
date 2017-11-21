//#include <MainLib.h> //change
#include "RandomGenerator.h" //change
#include <algorithm> 
#include <cassert> //change

#ifndef unix
#define CONSTANTS_VERY_LARGE_INT 4294967291
#define CONSTANTS_VERY_LARGE_INT_2 362436069
#define CONSTANTS_VERY_LARGE_INT_3 4294967295
#else
#define CONSTANTS_VERY_LARGE_INT  4294967291LL
#define CONSTANTS_VERY_LARGE_INT_2 362436069LL
#define CONSTANTS_VERY_LARGE_INT_3 4294967295LL
#endif

tRandomGenerator _RandomProbGenerator;

unsigned long
tMarsagliaGenerator::Next()
{
  register unsigned long hold;
  if (--i==0) i=43;
  if (--j==0) j=43;
  hold = word1[i]+carry;
  if (word1[j] < hold)
    {
      word1[i] = CONSTANTS_VERY_LARGE_INT - (hold-word1[j]);
      carry = 1;
    }
  else
    {
      word1[i] = word1[j]-hold;
      carry=0;
    }
  weyl-=CONSTANTS_VERY_LARGE_INT_2;
  return word1[i] - weyl;
}

tMarsagliaGenerator::tMarsagliaGenerator(unsigned query)
{
  unsigned start=1;  // default

  if (query)
    {
      cout << "Enter random number seed, 0<seed<65000." << endl;
      cin >> start;
    }

  srand(start);

  for (j=1;j<=43;j++)        		// make initial numbers.
    {
      word1[j]=rand();
      word1[j]=((word1[j] << 15 ) + rand());
      word1[j]=((word1[j] << 2 ) + rand()%4);
      if (word1[j]>CONSTANTS_VERY_LARGE_INT) word1[j] = CONSTANTS_VERY_LARGE_INT;
    }

  // initialize markers
  i=44;  j=23;  carry=1; weyl=rand();

  for (unsigned long a=1,garbage; a<100000; a++)       	// Warm-up
    garbage = Next();
}

void
tMarsagliaGenerator::Initialize(unsigned start)
{
  srand(start);

  for (j=1;j<=43;j++)        		// make initial numbers.
    {
      word1[j]=rand();
      word1[j]=((word1[j] << 15 ) + rand());
      word1[j]=((word1[j] << 2 ) + rand()%4);
      if (word1[j]>CONSTANTS_VERY_LARGE_INT) word1[j] = CONSTANTS_VERY_LARGE_INT;
    }

  i=44;  j=23;  carry=1; weyl=rand();				// initialize markers

  for (unsigned long a=1,garbage; a<100000; a++)       	// Warm-up
    garbage = Next();
}

double tMarsagliaGenerator::RandomDouble(double range)
{
  return ((((double)Next())/CONSTANTS_VERY_LARGE_INT_3)*range);
}



double
tRandomGenerator::SampleUniform(void)
{
  return RandomDouble(1.0);
}



double
tRandomGenerator::SampleNormal()
{
#ifdef notdef  
  // Ratio of uniforms, Ripley, Stochastic Simulation, p. 82
  double U, V;
  double X;
  double Z;
  double c = sqrt(2/exp(1));
  do {
  retry:
    U = c*(SampleUniform()*2-1);
    V = c*(SampleUniform()*2-1);
    if( fabs(U) < EPS )
      goto retry;
    X = V/U;
    Z = .25*X*X;
    if( Z < 1 - U )
      break;
  } while( (Z > 0.259/U + 0.35 ) || Z > -log(U) );
  return X;
#endif

  // Box-Muller, Ripley p. 54
  
  double theta = 2*3.1415926535897932384626433832795*SampleUniform(); //change
//  double theta = 2*M_PI*SampleUniform(); //change
  double R = sqrt(2*-log(SampleUniform()));
  return R * cos(theta);
}

// Code adopted from David Heckerman
//-----------------------------------------------------------
//	DblGammaGreaterThanOne(dblAlpha)
//
//	routine to generate a gamma random variable with unit scale and
//      alpha > 1
//	reference: Ripley, Stochastic Simulation, p.90 
//	Chang and Feast, Appl.Stat. (28) p.290
//-----------------------------------------------------------
double tRandomGenerator::DblGammaGreaterThanOne(double dblAlpha)
{
  double rgdbl[6];

  rgdbl[1] = dblAlpha - 1.0;
  rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
  rgdbl[3] = 2.0 / rgdbl[1];
  rgdbl[4] = rgdbl[3] + 2.0;
  rgdbl[5] = 1.0 / sqrt(dblAlpha);

  for (;;)
  {
    double  dblRand1;
    double  dblRand2;
    do
    {
      dblRand1 = SampleUniform();
      dblRand2 = SampleUniform();

      if (dblAlpha > 2.5)
	dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);

    } while (!(0.0 < dblRand1 && dblRand1 < 1.0));

    double dblTemp = rgdbl[2] * dblRand2 / dblRand1;

    if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
	rgdbl[3] * log(dblRand1) + dblTemp - log(dblTemp) < 1.0)
    {
      return dblTemp * rgdbl[1];
    }
  }
  assert(false);
  return 0.0;
} 


      
/* routine to generate a gamma random variable with unit scale and alpha
< 1

   reference: Ripley, Stochastic Simulation, p.88 */

double
tRandomGenerator::DblGammaLessThanOne(double dblAlpha)
{
  double dblTemp;

  const double	dblexp = exp(1.0);

  for (;;)
  {
    double dblRand0 = SampleUniform();
    double dblRand1 = SampleUniform();
    if (dblRand0 <= (dblexp / (dblAlpha + dblexp))) 
    {
      dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
		    dblexp, 1.0 / dblAlpha);
      if (dblRand1 <= exp(-1.0 * dblTemp)) 
	return dblTemp;
    }
    else 
    {
      dblTemp = -1.0 * log((dblAlpha + dblexp) * (1.0 - dblRand0) /
			   (dblAlpha * dblexp)); 
      if (dblRand1 <= pow(dblTemp,dblAlpha - 1.0)) 
	return dblTemp;
    }				
  }
  assert(false);
  return 0.0;
}  /* DblGammaLessThanOne */

// Routine to generate a gamma random variable with unit scale (beta = 1)
double
tRandomGenerator::SampleGamma(double dblAlpha)
{
  assert(dblAlpha > 0.0);
  if( dblAlpha < 1.0 )
    return DblGammaLessThanOne(dblAlpha);
  else
    if( dblAlpha > 1.0 )
      return DblGammaGreaterThanOne(dblAlpha);

  return -log(SampleUniform());
}  /* gamma */

//-------------------------------------------------------------------------
// SampleDirichlet expects an RGALPHA such that the Dirichlet's parameters
// are a vector of alphas such that
// p(\theta) = c * \prod_i  \theta_i^{\alpha_i - 1} 
//-------------------------------------------------------------------------

void
tRandomGenerator::SampleDirichlet(const vector<double>& mean,
				  double precision,
				  vector<double>& rgprob)
{
  assert(mean.size() == rgprob.size());
  double dblSum = 0.0;

  int i;
  
  for( i = 0; i < mean.size(); i++ )
    dblSum += (rgprob[i] = SampleGamma(mean[i] * precision));

  for( i = 0; i < rgprob.size(); i++)
  {
    rgprob[i] = rgprob[i] / dblSum;
  }

#ifdef LOG  
  cerr << "SampleDirichlet(";
  for(int j = 0; j < mean.size(); j++ )
    cerr << mean[j]*precision <<" ";
  cerr << ") = [";
  for(int j = 0; j < mean.size(); j++ )
    cerr << rgprob[j] <<" ";
  cerr <<"]\n";
#endif
} 

void
tRandomGenerator::SampleDirichlet(int n,
				  double const * mean,
				  double precision,
				  double * rgprob)
{
  double dblSum = 0.0;
  int i;
  
  for( i = 0; i < n; i++ )
    dblSum += (rgprob[i] = SampleGamma(mean[i] * precision));
  
  for( i = 0; i < n; i++)
  {
    rgprob[i] = rgprob[i] / dblSum;
  }

#ifdef LOG  
  cerr << "SampleDirichlet(";
  for(int j = 0; j < n; j++ )
    cerr << mean[j]*precision <<" ";
  cerr << ") = [";
  for(int j = 0; j < n; j++ )
    cerr << rgprob[j] <<" ";
  cerr <<"]\n";
#endif
} 

void
tRandomGenerator::SampleDirichlet(int n,
				  float const * mean,
				  double precision,
				  float * rgprob)
{
  double dblSum = 0.0;
  int i;
  
  for( i = 0; i < n; i++ )
    dblSum += (rgprob[i] = SampleGamma(mean[i] * precision));

  for( i = 0; i < n; i++)
  {
    rgprob[i] = rgprob[i] / dblSum;
  }
} 

int
tRandomGenerator::SampleMultinomial( int n, double const* prob )
{
  double p = SampleUniform();
  int i = 0;
  do{
    p -= prob[i++];
  } while( p > 0.0 && i < n );

#ifdef LOG  
  cerr << "SampleMulti(";
  for(int j = 0; j < n; j++ )
    cerr << prob[j] <<" ";
  cerr << ") = " << i-1 << "\n";
#endif
  
  return i-1;
}

int
tRandomGenerator::SampleLogMultinomial( int n, double* prob )
{
  int i;
  double s = 0.0;
  double M = prob[0];
  for( i = 1; i < n; i++ )
    if( prob[i] > M )
      M = prob[i];
  
  for( i = 0; i < n; i++ )
  {
    prob[i] = exp(prob[i] - M);
    s +=prob[i];
  }
  for( i = 0; i < n; i++ )
    prob[i] /= s;

  return SampleMultinomial(n , prob );
}



int
tRandomGenerator::sampleMulti( tMultinomial const & m ) {
  int n;double const * p = NULL;
  p = m.getParams(n);
  return SampleMultinomial(n,p);
}

vector<int>
tRandomGenerator::sampleGroup(int size,int groupSize)
{
  assert(groupSize<=size);

  vector<int> group;

  if (groupSize==size) {
    for( int i=0; i<groupSize; i++ )
      group.push_back(i);
    return group;
  }
  
  vector<int> v;
  for( int i=0; i<size; i++)
    v.push_back(i);
  
  for( int count = 0; count<groupSize; count++)
  {
    double val = SampleUniform();
    int chosenPos = (int)( val*(size-count) );
    group.push_back( v[chosenPos] );
    v[chosenPos] = v[size-1-count];
  }
  sort(group.begin(),group.end());

  return group;
}
