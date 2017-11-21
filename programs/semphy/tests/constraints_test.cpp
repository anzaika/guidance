#include "constraints.h"
#include "someUtil.h"

int main(int argc, char *argv[])
{
  if (argc<3) {
    cout << "usege: "<<argv[0]<<" input-tree-filename constraint-tree-filename"<<endl;
    exit (1);
  }
  cout<<"input tree"<<endl;
  string tf(argv[1]);
  tree t(tf);
  t.output(cout,tree::PHYLIP, true);
  cout<<endl<<endl;


  cout<<"constraint tree"<<endl;
  string cf(argv[2]);
  tree ct(cf);
  ct.output(cout);
  cout<<endl<<endl;
  
  constraints c(ct);
  c.setTree(t);
  if (c.fitsConstraints())
    cout<<"fits constraints"<<endl;
  else {
    cout<<"does not fit constraints:"<<endl;
    c.outputMissingClads(cout);
  }
  cout <<endl<<endl;

  vector<tree::nodeP> allNodes;
  t.getAllNodes(allNodes,t.getRoot());
  for (vector<tree::nodeP>::iterator i=allNodes.begin();i!=allNodes.end();++i)
    cout << (*i)->id()<< " "<<(*i)->name()<<endl;

cout <<endl<<endl;
  
  cout << "constraint table" <<endl;
  cout<<c.getPeneltyTable()<<endl;
  return (0);
}
