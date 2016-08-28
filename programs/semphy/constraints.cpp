// $Id: constraints.cpp 409 2005-06-28 13:12:24Z ninio $

#include "constraints.h"
#include "someUtil.h"

constraints::constraints(const string& constraingFileName):_ct(constraingFileName),_peneltyTable(_ct.getLeavesNum()*2-2)
{
  constraints_real_constractor();
};

constraints::constraints(const tree& in_ct):_ct(in_ct),_peneltyTable(_ct.getLeavesNum()*2-2)
{
  constraints_real_constractor();
}

void constraints::constraints_real_constractor()
{
  for(int i=0;i<_peneltyTable.size();++i)
    _peneltyTable[i].resize(_peneltyTable.size());
  tree::nodeP s=NULL;
  for( s=_ct.getRoot(); s->isInternal();s=s->getSon(0));
  _sonOfRoot=s->name();
  _ct.rootAt(s->father());	// move root to standard location
  make_list_of_constraints();
}

void test(constraints::csidp& a){};

constraints::csidp constraints::recursivlyFillTable(const tree::nodeP root)
{
  csidp csi;
  if (root->isLeaf()){		// end of recursion
    csi.first.insert(root->name());
    csi.second.insert(root->id());
  } else {
    for (int i=0;i < root->getNumberOfSons();++i){
      csidp csis(recursivlyFillTable(root->getSon(i)));

	{ for (childrenSet::iterator j=csis.first.begin(); j!=csis.first.end();++j) csi.first.insert(*j);}
	{ for (set<int>::iterator j=csis.second.begin(); j!=csis.second.end();++j) csi.second.insert(*j);}
      csi.second.insert(root->id());
    }
    childrenSetSet::iterator t=_cssTmp.find(csi.first);
    if (t != _cssTmp.end()){
      addPenetly(csi.second);
      _cssTmp.erase(csi.first);
    } 
  }
  return (csi);
}

void constraints::setTree(const tree& t1)	// give the current phlogony, recompute penerly table
{
  mult(_peneltyTable,0.0);
  _currentTree = t1;
  vector<tree::nodeP> allNodes;
  _currentTree.getAllNodes(allNodes,_currentTree.getRoot());
  for (vector<tree::nodeP>::iterator i=allNodes.begin();i!=allNodes.end();++i)
    if ((*i)->name()==_sonOfRoot) {
      _currentTree.rootAt((*i)->father());	// move root to standard location 
      break;
    }
  _cssTmp=_css;
  recursivlyFillTable(_currentTree.getRoot());
}

bool constraints::fitsConstraints()
{
  return(_cssTmp.empty());
}

void constraints::outputMissingClads(ostream& out)
{
  for (childrenSetSet::iterator i=_cssTmp.begin();i!=_cssTmp.end();++i){
    out <<"(";
    for (childrenSet::iterator j=i->begin();j!=i->end(); /*++j*/){
      out << *j;
      ++j;
      if (j != i->end())
	out <<",";
      else
	out <<")"<<endl;
    }
  }
}

void constraints::addPenetly(set<int>& v)
{
  for (set<int>::iterator i=v.begin();i!= v.end();++i){
    for(int k=0;k<_peneltyTable.size();++k){
      _peneltyTable[*i][k]++;
      _peneltyTable[k][*i]++;
    }
    for (set<int>::iterator j=v.begin();j!= v.end();++j){
      _peneltyTable[*i][*j]--;
      _peneltyTable[*j][*i]--;
    }
  }
}


const VVdouble& constraints::getPeneltyTable()
{
  return (_peneltyTable);
}


void constraints::make_list_of_constraints()
{
  recursive_make_list_of_constraints(_ct.getRoot());  
}

constraints::childrenSet constraints::recursive_make_list_of_constraints(const tree::nodeP& root)
{
  childrenSet cs;
  if (root->isLeaf()){		// end of recursion
    cs.insert(root->name());
  } else {
    for (int i=0;i<root->getNumberOfSons();++i){
      childrenSet scs(recursive_make_list_of_constraints(root->getSon(i)));
      for (childrenSet::iterator j=scs.begin(); j!=scs.end();++j) cs.insert(*j);
    }
    _css.insert(cs);
  }
  return (cs);
}

