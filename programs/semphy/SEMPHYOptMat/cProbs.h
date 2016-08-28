// cProbs.h: interface for the cProbs class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CPROBS_H__64290E26_D9E7_4687_AB56_6639348AD834__INCLUDED_)
#define AFX_CPROBS_H__64290E26_D9E7_4687_AB56_6639348AD834__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "cRow.h"

class cProbs : public cRow  
{
public:
	cProbs(const cRow& row);
	cProbs(size_t Size);	
  cProbs(istream& in);
	virtual ~cProbs(void);

};

#endif // !defined(AFX_CPROBS_H__64290E26_D9E7_4687_AB56_6639348AD834__INCLUDED_)
