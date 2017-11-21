#include "cSquareMatrix.h"
#include "cOccuranceData.h"
#include <string>
#include <stdio.h>

cOccurancePair ReadcOccurancePair(istream& in, const int length)
{
  cSquareMatrix m(length);
  double dist=0.0;
  string tmp;
  in >> tmp;
  sscanf(tmp.c_str(),"dist=%lf",&dist);
  m.Read(in);
  
  cOccurancePair p(dist, m);
  return(p);
}

