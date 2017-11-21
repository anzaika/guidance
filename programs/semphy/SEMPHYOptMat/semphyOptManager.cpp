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
#include "cOptManager.h"


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

  cOptManager::tOptParams params;

 //    params.OptType = cOptManager::PARALLEL;
    params.OptType = cOptManager::SEQUENTIAL;
  
  cOptManager Manager( D, Pi, params );

  cProbModel M = Manager();

  M.print(cout);
}

