#include "bootstrap.h"
#include "treeUtil.h"

#include "constantSemphyDistance.h"
#include "rearangeTree.h"



int main(int argc, char *argv[])
{
  cout << "creating a bootstrap object from a file"<<endl;
  
  string filename(argv[1]);

  vector<tree> tv(getStartingTreeVecFromFile(filename));

  // first constractor
  cout << " first constractor"<<endl;
  bootstrap b1(filename);
  b1.print();
  cout <<endl;

  cout << "geing whights from a tree" << endl;
  vector<float> v1(b1.GetWeightsToTree(tv[0])) ;
  for (vector<float>::iterator i = v1.begin();i!=v1.end();++i)
    cout << " "<<*i;
  cout << endl;

  cout << "checkout the support of a tree" <<endl;
  b1.printTreeNH(cout, tv[0], v1);
  cout <<endl <<endl;
  
  tree t(tv[0]);

    // add a "demi" distance
    semphyDistance* semDis1 = NULL;
    semDis1 = new constantSemphyDistance();
    

  // compatability tests
  {

    set< pair<int,int> > outset;
    vector<float> support;
    cout <<"compatability (0.5)"<<endl;
    support=b1.thresholdTree(t,outset);
    
    for (set< pair<int,int> >::iterator i=outset.begin();i!=outset.end();++i)
      {
	cout << "<"<<i->first<<","<<i->second <<">:"<<support[i->first]<<endl;
      }
    
    
    // note, you can not reuse the tree after the rearrange step!!!
    tree outtree1(t);		// copy
    
    b1.printTreeNH(cout, outtree1, support);
    cout<<endl;
    
    rearrangeTree  rt(&outset,semDis1);
    rt.reconstructTree(outtree1);
    
    outtree1.output(cout);
    b1.printTreeNH(cout, outtree1, support);
    cout<<endl;
    cout<<endl;
  }


  // compatability 0.1

{

    set< pair<int,int> > outset;
    vector<float> support;
    cout <<"compatability 0.1"<<endl;
    support=b1.thresholdTree(t,outset, 0.1);
    
    for (set< pair<int,int> >::iterator i=outset.begin();i!=outset.end();++i)
      {
	cout << "<"<<i->first<<","<<i->second <<">:"<<support[i->first]<<endl;
      }
    
    
    // note, you can not reuse the tree after the rearrange step!!!
    tree outtree2(t);		// copy
    
    b1.printTreeNH(cout, outtree2, support);
    cout<<endl;
    
    rearrangeTree  rt(&outset,semDis1);
    rt.reconstructTree(outtree2);
    
    outtree2.output(cout);
    b1.printTreeNH(cout, outtree2, support);
    cout<<endl;
  }
  // compatability 0.9

{

    set< pair<int,int> > outset;
    vector<float> support;
    cout <<"compatability 0.9"<<endl;
    support=b1.thresholdTree(t,outset, 0.9);
    
    for (set< pair<int,int> >::iterator i=outset.begin();i!=outset.end();++i)
      {
	cout << "<"<<i->first<<","<<i->second <<">:"<<support[i->first]<<endl;
      }
    
    
    // note, you can not reuse the tree after the rearrange step!!!
    tree outtree3(t);		// copy
    
    b1.printTreeNH(cout, outtree3, support);
    cout<<endl;
    
    rearrangeTree  rt(&outset,semDis1);
    rt.reconstructTree(outtree3);
    
    outtree3.output(cout);
    b1.printTreeNH(cout, outtree3, support);
    cout<<endl;
  }




  
//  if (semDis1!= NULL) delete semDis1;

  return (0);
}
