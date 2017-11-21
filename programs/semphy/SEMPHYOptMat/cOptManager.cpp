#include "StdAfx.h"
#include "cOptManager.h"
#include "cProbModel.h"
#include "cUpdateAngle.h"
#include "cEigenValsOptimise.h"


cOptManager::cOptManager( cOccuranceData const& Data, 
			  cProbs& BackgroundProbs,
			  tOptParams Params )
  : _params( Params ),
    _data( Data ),
    _initModel( BackgroundProbs )
{
}

cOptManager::cOptManager( cOccuranceData const& Data, 
			  cProbModel const& Model,
			  tOptParams Params )
  : _params( Params ),
    _data( Data ),
    _initModel( Model )
{
}

cOptManager::~cOptManager()
{
}

void
cOptManager::IterationAngleParallel(cDataProbModel& DM)
{
  cUpdateAngle Update(DM);
  cAngles NewA = Update.ComputeNewPhi();

  cout << "IterationAngleParallel" << endl;
  cout << "  Old = ";
  DM.ProbModel().Angles().print();
  
  // reset all frozen angles...
  for( size_t k = 0; k < First_Angle(); k++ )
    NewA.set( k, DM.ProbModel().Angles().get( k ) );

  DM.ProbModel().setAngles( NewA );
  
  cout << "  New = ";
  NewA.print();
  cout << "  New Angle LL = " << DM.LL() << endl;
}

double EmpiricalDerivAngle(cDataProbModel& DM, int k)
{
  double eps = 0.00000001;
  double origPhi=DM.ProbModel().Angles().get(k);
  DM.ProbModel().setAngle( k, origPhi-eps );
  double llMinus = DM.LL();
  DM.ProbModel().setAngle( k, origPhi+eps );
  double llPlus = DM.LL();
  DM.ProbModel().setAngle( k, origPhi );
  return (llPlus-llMinus)/(2*eps);
}

void
cOptManager::IterationAngleSequential(cDataProbModel& DM)
{
  cout << "IterationAngleSequential" << endl;

  // right now determinist round robin algorithm
  for( size_t k = First_Angle(); k < DM.ProbModel().Angles().NDims(); k++ )
  {
    double OldLL=DM.LL();
    double NewLL=OldLL;
    int IterNum=0;
    //    do {
    //      IterNum++;
      OldLL=NewLL;
      cUpdateAngle Update(DM);
      cAngles const& A  = DM.ProbModel().Angles();
      double phi = A.get(k);
      
      cout << "Emp Dphi "<<EmpiricalDerivAngle(DM,k)<<"\t";
      cout << "drotXY "  <<DM.DerivRotxy(A.GetDim(k).first, A.GetDim(k).second)<<"\t";
      
      double new_phi = Update.ComputeNewPhi(k);


      cout << "  Phi(" << k << ") ["
	   << A.GetDim(k).first << ", "
	   << A.GetDim(k).second << "] = "
	   << phi 
	   << " New = " << new_phi
	   << endl;
      
      DM.ProbModel().setAngle( k, new_phi );
      
      NewLL= DM.LL();
      
      if( isnan(NewLL) )
      {
	DM.ProbModel().setAngle( k, phi );
	NewLL = OldLL;
	cout << "Found nan! Reset to " << phi << endl;
	//	break;
      }
      if( NewLL < OldLL )
      {
	DM.ProbModel().setAngle( k, phi );
	cout << "Decrease in LL " << NewLL - OldLL << "! Reset to " << phi << endl;
	NewLL = OldLL;
	//	break;
      }

      cout << "  New Angle LL = " << NewLL << " \t Delta LL="<<NewLL-OldLL<< endl;
      //    }  while( isnan(NewLL)||isnan(OldLL)||(NewLL - OldLL > _params.StopCriteria*50 &&
      //					   IterNum < _params.MaxIter/10 ));

  } 
}

void
cOptManager::IterationEigen(cDataProbModel& DM)
{
  for (int x=1;x<DM.ProbModel().size();++x) {
    cEigenValsOptimise EVO(DM);

    cout << "  Eig("<<x<<") = "<< DM.ProbModel().GetD(x);

    double new_dx=EVO.optimize(x,0.001);

    cout << "  New = "<<new_dx <<endl;

    DM.ProbModel().SetD(x,new_dx);

    cout << "  New Eigen LL = " << DM.LL() << endl;
  }
}

cProbModel
cOptManager::operator()(void)
{
  cProbModel M( _initModel );
  cDataProbModel DM( _data, M );

  double OldLL;
  double NewLL = DM.LL();
  int IterNum = 0;
  
  do
  {
    IterNum++;
    OldLL = NewLL;

    switch( _params.OptType )
    {
    case PARALLEL:
      IterationAngleParallel(DM);
      break;
      
    case SEQUENTIAL:
      IterationAngleSequential(DM);
      break;
    };

    NewLL = DM.LL();
    cout << "AngleOpt[ " << IterNum << " ] = " << NewLL << endl;
    
    IterationEigen(DM);

    NewLL = DM.LL();
    cout << "EigenOpt[ " << IterNum << " ] = " << NewLL << endl;
    
  } while( NewLL - OldLL > _params.StopCriteria &&
	   IterNum < _params.MaxIter );

  return M;
}

size_t cOptManager::First_Angle() const
{
  if( _params.FreezeBackground )
    return _initModel.size()-1;
  else
    return 0;
}

