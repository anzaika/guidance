// SEMPHYOptMat.cpp : Defines the entry point for the console application.
//

#include "StdAfx.h"
#include "cOccuranceData.h"
#include "cOrthonormalMatrix.h"
#include <fstream>
#include "cAngles.h"
#include "cCoeffMatrices.h"
#include "cDataProbModel.h"
#include "cUpdateAngle.h"
#include "cEigenValsOptimise.h"

cSquareMatrix
EmpiricalDuax( cProbModel& M, cOccuranceData const& D )
{
  double eps = 0.0000000001;
  cDataProbModel DM(D, M);
  double LL = DM.LL();
  cSquareMatrix Res( M.size(), 0.0 );
  for( size_t a = 0; a < M.size(); a++ )
    for( size_t b = 0; b < M.size(); b++ )
    {
      cSquareMatrix U( M.GetU() );

      U.Set(a,b,U[a][b]+eps);
      double LL1 = DM.LL(U);
      Res.Set(a,b,(LL1-LL)/eps);
    }

  return Res;
}

cSquareMatrix
EmpiricalDuaxNonConst( cProbModel& M, cOccuranceData const& D )
{
  double eps = 0.0000000001;
  cDataProbModel DM(D, M);
  double LL = DM.NonConstLL();
  cSquareMatrix Res( M.size(), 0.0 );
  for( size_t a = 0; a < M.size(); a++ )
    for( size_t b = 0; b < M.size(); b++ )
    {
      cSquareMatrix U( M.GetU() );

      U.Set(a,b,U[a][b]+eps);
      double LL1 = DM.NonConstLL(U);
      Res.Set(a,b,(LL1-LL)/eps);
    }

  return Res;
}

cRow
EmpiricalDx( cProbModel& M, cOccuranceData const& D )
{
  double eps = 0.000001;
  cDataProbModel DM(D, M);
  double LL = DM.LL();
  cRow Res( M.size() );
  
  for( size_t a = 0; a < M.size(); a++ )
  {
    cProbModel M1 = M;

    M1.SetD( a, M.GetD(a) - eps );
    cDataProbModel DM1(D, M1);
    double LL1 = DM1.LL();
    Res[a] = (LL-LL1)/eps;
    
  }
  return Res;
}

cRow
EmpiricalDerivAngles(cProbModel& M, cOccuranceData const& D )
{
  double eps = 0.000001;
  cAngles const& A = M.Angles();
  cRow Res( A.NDims() );
  cDataProbModel DM(D, M);
  double LL = DM.LL();

  for( size_t i = 0; i < A.NDims(); i++ )
  {
    cProbModel M1 = M;
    M1.setAngle(i, A.get(i)-eps );
    cDataProbModel DM1(D, M1);
    double LL1 = DM1.LL();
    Res[i] = (LL-LL1)/eps;
    
  }
  return Res;
}

cSquareMatrix BuildMUT( cProbModel& M, double t )
{
  cSquareMatrix A( M.size(), 0.0 );
  for( size_t a = 0; a < M.size(); a++ )
    for( size_t b = 0; b < M.size(); b++ )
      A.Set(a,b, M.Muabt( a, b, t ));

  return A;
}

int main(int argc, char* argv[])
{
  if (argc<=2) 
  {
    cerr <<"usage: "<<argv[0]<<" data-file prob-file"<<endl;
    exit(1);
  }
  ifstream in, in2;
  in.open(argv[1]);
  
  cOccuranceData D(in);
  in.close();


  in2.open(argv[2]);
  cProbs Pi(in2);
  in2.close();

  cout<<"Read Pi: ";
  Pi.print();

  cProbModel M(Pi);

  cout<<"first Pi:" ;
  M.GetPi().print();  
  M.Angles().print();


  // Set diagonal
//   M.SetD(1,-0.3620);
//   M.SetD(2,-0.5078);
//   M.SetD(3,-0.2316);

  cout<<endl<<"first D:  " ;
  M.GetD().Get().print();

  cDataProbModel DM(D, M);


  cout<< "U"<<endl;
  M.GetU().print();

  cout<<"Q [Mut()]"<<endl;
  M.Mut(1.0).print();

  cout<<"Q [Muabt]"<<endl;
  BuildMUT(M,1.0).print();

  cout << "LL = " << DM.LL() << endl;
  cout << "NonConstLL = " << DM.NonConstLL() << endl;

  cout<< "direct duxy"<<endl;
  DM.DerivAllUxy().print();

  cout<< "Emp duxy"<<endl;
  EmpiricalDuax( M, D ).print();

  cout<< "Emp NonConst duxy"<<endl;
  EmpiricalDuaxNonConst( M, D ).print();

  cout<< "direct drphi"<<endl;
  DM.DerivAllAngles().print();

  cout<< "emp drphi"<<endl;
  EmpiricalDerivAngles( M, D ).print();
  
  cout<< "direct dDx"<<endl;
  DM.DerivAllEigen().print();
  
  cout<< "emp dDx"<<endl;
  EmpiricalDx( M, D ).print();
  
  cout<< "U1*U1'^t"<<endl;
  (M.GetU()*M.GetU().Transpose()).print();

  cout<<"Q"<<endl;
  M.GetQ().print();

  cout<<"R"<<endl;
  M.GetR().print();

  //  M.setAngle(5, 1.0);
  cAngles A(M.Angles());
  A.print();

  double oldLL=DM.LL();

  // Debug code
  for( int k = 0; k < A.NDims(); k++ )
  {
    cout << "Coeff(" << k << ") [" << A.GetDim(k).first << ", "
	 << A.GetDim(k).second << "]" << endl;
    
    cCoeffMatrices CoeffMatrices( A, k );
    cSquareMatrix M1(CoeffMatrices.CoeffMatrix(cDRotationMatrix::ONE));
    cSquareMatrix MS(CoeffMatrices.CoeffMatrix(cDRotationMatrix::SIN));
    cSquareMatrix MC(CoeffMatrices.CoeffMatrix(cDRotationMatrix::COS));
    
    double phi = A.get(k);
    for( size_t a = 0; a < M.size(); a++ )
      for( size_t b = 0; b < M.size(); b++ )
	M1.Set(a,b, M1[a][b]+ sin(phi)*MS[a][b] + cos(phi)*MC[a][b]);

    cout << "Reconstructed M (should be = U)\n";
    M1.print();
  }
  
  for (int i=0;i<40;++i) {
    //    cout<<"Mut()"<<endl;
    //    M.Mut(1.0).print();
    //    int k=(i/4)%(A.NDims()-A.size()+1)+A.size()-1;
    //    cout << "i="<<i<<" start nonconst = "<<DM.LL()<<"\t"<<DM.DerivRotxy(2,3)<<endl;
    //          k=5;

    cUpdateAngle UA(DM);	// need a new one!! derivatives are computed only once.

    cAngles newA=UA.ComputeNewPhi();		// compute all
    //    cAngles newA=UA.ComputeNewPhiNoPiChange();    // don't conpute fixed

    //    cout << "Opt(" << k << ") [" << A.GetDim(k).first << ", "
    //	 << A.GetDim(k).second << "]" << endl;

    //    double new_phi_k=UA.ComputeNewPhi(k);		// compute just one
    cout<<"new ";
    newA.print();
    M.setAngles(newA);
    //    double oldPhi=M.Angles().get(k);

    //    cout <<"phi("<<k<<") = "<<new_phi_k<<"\t";

	//    M.setAngle(k, new_phi_k);
    //     cout<< "U2"<<endl;
    //    
    double LL=DM.LL();
    cout << "i="<<i<<" end   LL = "<<LL<<"\t"<<LL-oldLL<<"\t"<< endl;
    M.Angles().print();
    cout<< "emp drphi"<<endl;
    EmpiricalDerivAngles( M, D ).print();

    oldLL=LL;
      
    cout << "U"<<endl;
    M.GetU().print();
    cout<<endl;
    
     //
//     cout<< "U2*U2^t"<<endl;
//     (M.GetU()*M.GetU().Transpose()).print();
    
//     cout<<"Q"<<endl;
//     M.GetQ().print();

//    M.GetPi().print();

    //    cout<<"R"<<endl;
    //   M.GetR().print();
  }

  cout<<endl<<" D:  " ;
  M.GetD().Get().print();


  cout << "U2"<<endl;
  M.GetU().print();

  cout<< "U2*U2^t"<<endl;
  (M.GetU()*M.GetU().Transpose()).print();
    
  cout<<"Q"<<endl;
  M.GetQ().print();

  cout <<"pi ";
  M.GetPi().print();

  cout<<"R"<<endl;
  M.GetR().print();

  cout <<"        LL:  "<<DM.LL()<<endl;
  cout <<"NonConstLL:  "<<DM.NonConstLL()<<endl;

  cout << endl<<"- - DIAGONAL - -"<<endl;
  for (int i=0;i<20;++i) {
    cEigenValsOptimise EVO(DM);
        //    int k=(i/4)%(A.NDims()-A.size()+1)+A.size()-1;
    int x=1+i%(A.size()-1);
    cRow EIG(DM.ProbModel().GetD().Get());
 
    
    cout <<"EIG: ";
    EIG.print();
    double dx=EIG[x];
    cout <<"old d("<<x<<")="<<dx<<endl;
//     for ( double tmp_dx = dx - 0.1 ; tmp_dx < dx + 0.1 ; tmp_dx += 0.001 ) {
//       cProbModel M2(M);
//       M2.SetD(x, tmp_dx);
//       cDataProbModel DM2(D,M2);
//       //cout << "derivative DerivEigen" << DM2.DerivEigen(x) << endl;      
//       //      cRow EIG2(DM2.ProbModel().GetD().Get());
//       //      cout << "tmpx " << tmp_dx << " ddx " << EVO.dLL_ddx(EIG2, x)<< " LL " << EVO.LL_dx(EIG2) << "\n";
//    }

    double lldx= EVO.LL_dx(EIG);
    double dlldx= EVO.dLL_ddx(EIG, x);
    cout <<"old lldx="<<lldx <<"  dlldx="<<dlldx<<endl;
    double new_dx=EVO.optimize(x,0.001);
    EIG[x]=new_dx;
    double nlldx=     EVO.LL_dx(EIG);
    
    cout <<"new lldx="<<nlldx<<"\t delta="<<nlldx-lldx<<"\tDllDdx = "<<EVO.dLL_ddx(EIG, x)<<endl;
    M.SetD(x,new_dx);
  }

  cout << "U2"<<endl;
  M.GetU().print();

  cout<< "U2*U2^t"<<endl;
  (M.GetU()*M.GetU().Transpose()).print();
    
  cout<<"Q"<<endl;
  M.GetQ().print();

  M.GetPi().print();

  cout<<"R"<<endl;
  M.GetR().print();


  return 0;
}

