#include "cEigenValsOptimise.h"
#include "numRec.h"
#include <cmath>

// double cEigenValsOptimise::ft_x(const int Which_t, const int x) const {
//   const cSquareMatrix& Mut=_DM.ProbModel().Mut(_DM.OccuranceData().Time(Which_t));
//   const cSquareMatrix& U=_DM.ProbModel().GetU();
//   const cSquareMatrix& S=_DM.OccuranceData().Matrix(Which_t);

//   double f=0.0;
  
//   for (int a=0;a<_DM.OccuranceData().alphabetSize();++a)
//     for (int b=0;b<_DM.OccuranceData().alphabetSize();++b) 
//       f += S[a][b] * U[x][a]*U[x][b]
// 	             /  Mut[a][b];
//   return f;
// }

// cRow cEigenValsOptimise::f_x(const int x) const {
//   cRow fx(_DM.OccuranceData().size(),0.0);
//   for (int Which_t=0;Which_t<_DM.OccuranceData().size();++Which_t)
//     fx[Which_t]=ft_x(Which_t,x);
//   return fx;	double DerivEigen(size_t i) const;
// }

double cEigenValsOptimise::LL_dx(const cRow& D) const
{
//   cout <<"current D ";
//   D.print(cout);
  double sum=0.0;
  for (int Which_t=0;Which_t<_DM.OccuranceData().size();++Which_t){
    const cSquareMatrix& S=_DM.OccuranceData().Matrix(Which_t);
    const cSquareMatrix Mutd=_DM.ProbModel().MutD(_DM.OccuranceData().Time(Which_t), D);
    //    Mutd.print(cout);
    for (int a=0;a<_DM.OccuranceData().alphabetSize();++a)
      for (int b=0;b<_DM.OccuranceData().alphabetSize();++b) 
	sum += S[a][b]*log(Mutd[a][b]);
  }
  return sum;
}

double cEigenValsOptimise::LL_dx(const double& dx, size_t x) const
{
  cRow D(_DM.ProbModel().GetD().Get());
  D[x]=dx;
  return(LL_dx(D));
}

double cEigenValsOptimise::dLL_ddx(const cRow& D, const int x) const
{
  const double dx=D[x];
  const cSquareMatrix& U=_DM.ProbModel().GetU();

  double sum=0.0;
  for (int Which_t=0;Which_t<_DM.OccuranceData().size();++Which_t){
    double t=_DM.OccuranceData().Time(Which_t);
    const cSquareMatrix& S=_DM.OccuranceData().Matrix(Which_t);
    const cSquareMatrix Mutd=_DM.ProbModel().MutD(t, D);
    double t_sum=0.0;
    
    for (int a=0;a<_DM.OccuranceData().alphabetSize();++a) 
      for (int b=0;b<_DM.OccuranceData().alphabetSize();++b)
	t_sum += S[a][b]*U[x][a]*U[x][b]/(Mutd[a][b]);
    //	cout << "EVO : a " << a << " b " << b << " add " << t*exp(t*dx)*S[a][b]*U[x][a]*U[x][b]/(Mutd[a][b]) << " addshort " << S[a][b]*U[x][a]*U[x][b]/(Mutd[a][b]) << " mut " << Mutd[a][b] << " sum " << sum + t*exp(t*dx)*t_sum << endl;
    //}
    //}
    sum += t*exp(t*dx)*t_sum;
  }
  
  return sum;
}


class C_EigenVals{
private:
  const cEigenValsOptimise& _cEO;
  cRow _d;
  const int _x;
public:
  C_EigenVals(const cEigenValsOptimise& cEO, const cRow& d, const int x):_cEO(cEO),_d(d),_x(x){};

  double operator() (double dx) {
    _d[_x]=dx;
    return(_cEO.LL_dx(_d));
  }
};

class C_EigenVals_d{
private:
  const cEigenValsOptimise& _cEO;
  cRow _d;
  const int _x;
public:
  C_EigenVals_d(const cEigenValsOptimise& cEO, const cRow& d, const int x):_cEO(cEO),_d(d),_x(x){};

  double operator() (double dx) {
    _d[_x]=dx;
    return(_cEO.dLL_ddx(_d,_x));
  }
};

template <typename regF, typename dF>
double _optimize_gradient(double old_x, regF f, dF df, 
					      const double init_delta, const double delta_factor,
					      const double min_step) 
{
  double step = init_delta;
  double new_x = old_x;
  double old_f, new_f;
  do {
    old_x = new_x;
    old_f = f(old_x);
    new_x = _improve_gradient(old_x, f, df, step);
    new_f = f(new_x);
    step *= delta_factor;
  } while ( step > min_step );
  return new_x;
}

template <typename regF, typename dF>
double _improve_gradient(double old_x, regF f, dF df, const double delta)
{
  double new_x = old_x;
  double new_f, old_f;
  do {
    old_x = new_x;
    old_f = f(old_x);
    //    double ddf = (df(old_x + delta) - df(old_x - delta))/2*delta;
    double dfoldx(df(old_x));
    if ( dfoldx > 0 ) {
      new_x = old_x + delta;      
    } else {
      new_x = old_x - delta;
    }
    //    new_x = old_x - df(old_x)/ddf;
    new_f = f(new_x);
  } while ( new_f > old_f + delta );
  return old_x;
}


double cEigenValsOptimise::optimize(const int x, const double tol) const 
{
  double newdx=_optimize_gradient( _DM.ProbModel().GetD(x), 
				  C_EigenVals(*this,_DM.ProbModel().GetD().Get(),x),
				  C_EigenVals_d(*this,_DM.ProbModel().GetD().Get(),x), .01, .1, .000001 );
				  
   //   double resL = -dbrent(MAXEIGENVAL,_DM.ProbModel().GetD(x),MINEIGENVAL,
   //		 C_EigenVals(*this,_DM.ProbModel().GetD().Get(),x),
   //		 C_EigenVals_d(*this,_DM.ProbModel().GetD().Get(),x),
   //		 tol,
   //		 &dist);
  double resL = LL_dx(newdx, x);
   cout <<" after hillclimb LLest="<<resL<<endl;
   return newdx;
}
