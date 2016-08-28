// cOccuranceData.h: interface for the cOccuranceData class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COCCURANCEDATA_H__B97B0121_A0A0_4E6C_A9E3_B0914CA46423__INCLUDED_)
#define AFX_COCCURANCEDATA_H__B97B0121_A0A0_4E6C_A9E3_B0914CA46423__INCLUDED_

#include "cSquareMatrix.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <vector>
using namespace std;
#include "cOccurancePair.h"

typedef vector<cOccurancePair> tOccuranceData;

class cOccuranceData
{
public:
	size_t alphabetSize(void) const;
	void AddPair(const cOccurancePair& NextPair);
	cSquareMatrix AverageNormalization(void);
	cOccuranceData(const cOccurancePair& FirstPair);
	cOccuranceData(const int size);
	cOccuranceData(istream& in);
	virtual ~cOccuranceData();
	void print(ostream& sout = cout) const;
  const cSquareMatrix& Matrix(const int i) const {return _data[i].Matrix();};
  
  
  // return S^(timeid)_[a,b] + S^(timeid)_[b,a]
  inline double getSS(const int time_id,const int a,const int b) const;
  inline const double& getSstar(const int a) const;
  inline double Time(const int time_id) const {return _data[time_id].Time();};
  inline size_t size() const {return _data.size();};
  //  double SstartA(const int a);
  //  double SstartB(const int b);
private:
  tOccuranceData _data;
	size_t _mat_size;
  vector<double> _SstarDelta;	// S^*_[*,y] - S^*_[y,*]  - IS THIS RIGHT???
  //  vector<vector<vector<double> > > _StabSUM;	// S^t_[a,b] + S^*_[b,a]

  void _updateSuffisiontStatistics(void);

};

#endif // !defined(AFX_COCCURANCEDATA_H__B97B0121_A0A0_4E6C_A9E3_B0914CA46423__INCLUDED_)
