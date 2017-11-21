// 	$Id: Multinomial.cpp 414 2005-06-29 16:40:37Z ninio $	


//#include <MainLib.h> 
#include "Multinomial.h"// CHANGE
#include <vector>
#include <cmath> // CHANGE
#include <cassert> // CHANGE
using namespace std;// CHANGE

//int tMultinomial::MaxN = 4;

// void
// tMultinomial::setMaxN(int n) {
//   assert(n>0);
//   MaxN = n;
// }

tMultinomial::tMultinomial()
{
  _n = 0;
}

tMultinomial::tMultinomial(int n)
{
  _n = n;
  for( int i = 0; i < n; i++ )
    _PR1[i] = 0;
}

tMultinomial::tMultinomial(const tMultinomial& P )
{
  _n = P._n;
  for( int i = 0; i < _n; i++ )
    _PR1[i] = P._PR1[i];

}

tMultinomial::tMultinomial(const vector<double>& P,const bool isExp )
{
  _n = P.size();
  for( int i = 0; i < _n; i++ )
    _PR1[i] = P[i];
  isExp ? NormalizeExp() : Normalize();
}

tMultinomial
tMultinomial::operator=(const tMultinomial& P )
{
  _n = P._n;
  for( int i = 0; i < _n; i++ )
    _PR1[i] = P._PR1[i];
  return *this;
}


void
tMultinomial::Normalize( )
{
  int i;
  double T = 0;
  for(  i = 0; i < _n; i++ )
    T +=_PR1[i];
  
  if( T < 0.00001 )
    for(  i = 0; i < _n; i++ )
      _PR1[i] = 1.0/_n;
  else
    for(  i = 0; i < _n; i++ )
      _PR1[i] /= T;
}

void
tMultinomial::NormalizeExp( )
{
  int i;

  double M = _PR1[0];
  for(  i = 1; i < _n; i++ )
    if( _PR1[i] > M )
      M = _PR1[i];

  double T = 0;
  
  for(  i = 0; i < _n; i++ )
  {
    _PR1[i] = exp( _PR1[i] - M );
    T +=_PR1[i];
  }
  
  if( T < 0.00001 )
    for(  i = 0; i < _n; i++ )
      _PR1[i] = 1.0/_n;
  else
    for(  i = 0; i < _n; i++ )
      _PR1[i] /= T;
}


void tMultinomial::NormalizeCompExp()
{
  NormalizeExp();
  for (int  i = 0 ; i < _n; i++ )
    _PR1[i] = 1 - _PR1[i];
  Normalize();
}

/**
 * its the user reponsibility to use it on normallized dist only.
 * also note the user should avoid using it with prob of 0 in one of the values.
 * the multinomials themselves should be of the same size. 
 */

double
tMultinomial::KL(tMultinomial const & m ) const {
  assert (m._n == _n);
  double val = 0;
  for( int i = 0; i < _n; i++ )
  {
    val += (_PR1[i]* (log(_PR1[i])-log(m._PR1[i])));

  }
  return val;
}

/* CHANGE - NOT COMPILING IN VC++
ostream& operator<<(ostream & o, tMultinomial const& P)
{
  o << "M[ ";
  for( int i =0 ; i < P._n; i++ )
  {
    if( i > 0 )
      o << " ";
    o << P._PR1[i];
  }
  o << " ]";
  return o;
}
*/
void
tMultinomial::Read(istream& in )
{
  string s;
  in>>s;
  assert(!s.compare("M["));
  for( int i = 0 ; i < _n; i++ )
  {
    in>>_PR1[i];
  }
  in>>s;
  assert(!s.compare("]"));
  /*
  char *t;
  t = NextToken(in);
  assert( !strcmp(t, "M[" ) );
  for( int i = 0 ; i < _n; i++ )
  {
    char *t = NextToken(in);
    _PR1[i] = atof(t);
    //    cerr << t << " -> " << _PR1[i] << "\n";
  }
  
  t = NextToken(in);
  assert( !strcmp(t, "]" ) );
  */
}

int
tMultinomial::MAP() const
{
  int M = 0;
  for( int i = 1 ; i < _n; i++ )
    if( _PR1[i] > _PR1[M] )
      M = i;

  return M;
}
