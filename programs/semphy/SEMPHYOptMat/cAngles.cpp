#include "cAngles.h"
#include <cmath>

cAngles::cAngles(size_t size)
: _size(size), _phi(0), _dim(0)
{
	_init_phi_dim();
}



cRotationMatrix cAngles::ROTk(const size_t k) const
{
	tPlane CurrentPlane(GetDim(k));
	return cRotationMatrix(_size, CurrentPlane.first,
		CurrentPlane.second, _phi[k]);
}

void cAngles::RotateMatrix(cOrthonormalMatrix& mat) const
{
  // multiply last first
  RotateMatrix(mat, _lastK(),-1);

    // MATAN: I don't like this,
    // but one has to do this if
    // one wants to get plain "0" to rotate!!!

}

size_t cAngles::NDims() const
{
	return _size*(_size-1)/2;
}

size_t cAngles::_lastK() const
{
	return NDims()-1;
}

void cAngles::RotateMatrix(cOrthonormalMatrix &mat, int kfirst, int klast) const
{
  //  cout << "Rotating" << endl;
  cOrthonormalMatrix tmp(mat);
  if (kfirst==klast) return;
  int step = -1;
  if ( kfirst < klast )
    step = 1;
  for( size_t k = kfirst ; k != klast ; k += step ) {			
    tPlane CurrentPlane(GetDim(k));
    mat.RotateAlpha(CurrentPlane.first, CurrentPlane.second, _phi[k]);

//  		cout << "Rot("<<k<<")=<"<<CurrentPlane.first<<","<<CurrentPlane.second<<"> by "<< _phi[k]<<endl;
// 		cout << "mat"<< endl;
// 		mat.print(cout);
// 		cout << "mat*mat'" << endl;
//  		(mat*mat.Transpose()).print(cout);

//  		tmp.Unit();
//  		cout<<"tmp"<<endl;
//  		tmp.print();
//  		tmp.RotateAlpha(CurrentPlane.first, CurrentPlane.second, _phi[k]);
//  		cout<<"tmp.rot"<<endl;
//  		tmp.print();
//  		cout<<"u*u^t"<<endl;
//  		(tmp*tmp.Transpose()).print(cout);		
  }
}

size_t cAngles::size() const
{
	return _size;
}

void cAngles::print(ostream& sout) const
{
  sout<<"Angles: ";
  for (int i=0;i<NDims();++i){
    sout << _phi[i];
    sout <<" ";
  }
  sout << endl;
}

cAngles::tPlane cAngles::GetDim(size_t k) const
{
	return _dim[k];
}

cAngles::cAngles(const cProbs &BackgroundProbs)
  : _size(BackgroundProbs.size()), _phi(0), _dim(0)
{
  // Debugged by Nir & Matan
  
  _init_phi_dim();
  double factor = 1.0;
  for ( size_t i = size()-1 ; i >0 ; --i )
  {
    _phi[i-1] = -asin(-sqrt(BackgroundProbs[i])/factor);
#ifdef debugcode
    cout << "i: "<< i << " _phi " << _phi[i-1] ;
    cout << " pi " << BackgroundProbs[i] << " " << sqrt(BackgroundProbs[i]) << endl;
    cOrthonormalMatrix M(size());
    M.Unit();
    RotateMatrix(M);
    M.print();
    cout << endl; 
#endif
    factor *= cos( _phi[i-1] );
    //    cout << "Factor = " << factor << endl;
  }
}

void cAngles::_init_phi_dim()
{
	_dim.reserve(NDims());
	_phi.reserve(NDims());
	for ( size_t i = 0 ; i < _size-1 ; ++i ) {
		for ( size_t  j = i+1 ; j < _size ; ++j ) {
			tPlane CurrentPlane(i,j);
			_dim.push_back(CurrentPlane);
			_phi.push_back(0.);
		}
	}
// 	for (int k=0;k<NDims();++k)
// 	  cout << "dim("<<k<<")= ["<<_dim[k].first<<"/"<<_dim[k].second<<"]"
// 	       <<endl;

}

size_t
cAngles::AngleIndex(int i, int j ) const
{
  int k; 
#ifdef notdef  
  k = j-i-1;
  for( int a = i-1 ; a >= 0; a-- )
    k += size()-a-1;

  cout << "i = " << i;
  cout << " j = " << j;
  cout << " k = " << k;
  cout << " k' = " << j-1+i*(size()-1)-(i*(i+1))/2;
  cout << endl;
#endif

  k = j-1+i*(size()-1)-(i*(i+1))/2;

  assert( GetDim(k).first == i );
  assert( GetDim(k).second == j );
  
  return k;
}

cAngles::cAngles( cOrthonormalMatrix const &U )
  : _size(U.size()), _phi(0), _dim(0)
{
  _init_phi_dim();

  cSquareMatrix Target(U);
  
  double factor = 1.0;
  for( size_t j = 0; j < size()-1 ; j++ )
  {
    cOrthonormalMatrix L(size());
    L.Unit();
    RotateMatrix(L);
    Target = L.Transpose()*U;
    factor = 1.0;
    
    for( size_t i = size()-1; i > j; i-- )
    {
      size_t k = AngleIndex( j, i );
      _phi[k] = -asin(Target[i][j]/factor);
      factor *= cos( _phi[k] );

#ifdef notdef
      cout << "k = " << k;
      cout << " j,i: "<< j << ", " << i << " _phi " << _phi[k] ;
      cout << " U[i,j] " << U[i][j] << endl;
      
      cOrthonormalMatrix M(size());
      M.Unit();
      RotateMatrix(M);
      M.print();
      cout << endl; 
#endif
    }
  }
}
