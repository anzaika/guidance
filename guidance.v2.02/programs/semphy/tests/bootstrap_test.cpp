#include "bootstrap.h"
#include "treeUtil.h"

#include "constantSemphyDistance.h"
#include "rearangeTree.h"

typedef map<int,MDOUBLE> Wmap;


int main(int argc, char *argv[])
{
  cout << "creating a bootstrap object from a file"<<endl;
  
  string filename("bootstrap_test.txt");
  if (argc>1)
    filename=argv[1];
  vector<tree> tv(getStartingTreeVecFromFile(filename));

  // first constractor
  cout << " first constractor"<<endl;
  bootstrap b1(filename);
  b1.print();
  cout <<endl;

  // secound constractor
  cout << " secound constractor" <<endl;
  bootstrap b2(tv);
  b2.print();
  cout <<endl;
  
  cout << "geing whights from a tree" << endl;
  Wmap v1(b1.getWeightsForTree(tv[0])) ;
  for (Wmap::iterator i = v1.begin();i!=v1.end();++i)
    cout << " "<<i->second;
  cout << endl;

  cout << "checkout the support of a tree" <<endl;
  b1.printTreeWithBPvalues(cout, tv[0], v1, false);
  cout <<endl <<endl;
  

  cout<< "remove the first tree from the list, and use is as bases for additional computation"<<endl;
  tree t(tv[0]);
  tv[0]=tv[1];
  
  // secound constractor
  bootstrap b3(tv);
  b3.print();
  
  Wmap  v3(b3.getWeightsForTree(t)) ;
  for (Wmap::iterator i = v3.begin();i!=v3.end();++i)
    cout << " "<<i->second;
  cout << endl;
  
  cout << "checkout the support of the removed tree"<<endl;
  b3.printTreeWithBPvalues(cout, t, v3, false);
  cout <<endl <<endl;


    // add a "demi" distance
    semphyDistance* semDis1 = NULL;
    semDis1 = new constantSemphyDistance();
    


  // compatability tests
  {
    cout <<"compatability 0.5 = Majority Rule"<<endl;
    tree outtree1(b1.consensusTree());
    cout << endl<<"from the consensus" <<endl;
    outtree1.output(cout);
     cout << "with support" <<endl;
     Wmap support(b1.getWeightsForTree(outtree1));
     
     b1.printTreeWithBPvalues(cout, outtree1, support, false);
    cout<<endl;
    cout<<endl;
  }



  // compatability 0.1
{
    cout <<"compatability 0.1"<<endl;
    tree outtree2(b1.consensusTree(0.1));
    cout << endl<<"from the consensus" <<endl;
    outtree2.output(cout);
     cout << "with support" <<endl;
     Wmap support(b1.getWeightsForTree(outtree2));
     
     b1.printTreeWithBPvalues(cout, outtree2, support, false);
    cout<<endl;
    cout<<endl;
}
  // compatability 0.9
{
    cout <<"compatability 0.9"<<endl;
    tree outtree3(b1.consensusTree(0.9));
    cout << endl<<"from the consensus" <<endl;
    outtree3.output(cout);
     cout << "with support" <<endl;
     Wmap support(b1.getWeightsForTree(outtree3));
     
     b1.printTreeWithBPvalues(cout, outtree3, support, false);
    cout<<endl;
    cout<<endl;
}
  
  // compatability 1.0
{
    cout <<"strict consensus"<<endl;
    tree outtree4(b1.consensusTree(1.0));
    cout << endl<<"from the consensus" <<endl;
    outtree4.output(cout);
     cout << "with support" <<endl;
     Wmap support(b1.getWeightsForTree(outtree4));
     
     b1.printTreeWithBPvalues(cout, outtree4, support, false);
    cout<<endl;
    cout<<endl;
}


  
  if (semDis1!= NULL) delete semDis1;

  return (0);
}
