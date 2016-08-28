#include "tree.h"
#include "treeUtil.h"

#include "someUtil.h"
#include "treeUtil.h"
#include "phylipFormat.h"
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;



void removeTaxa(tree& inTree,  const string& missingListFile);


int main(int argc, char **argv) {
	//talRandom::setSeed(2);
	
	if (argc != 4)
		errorMsg::reportError("input arguements should be: [1] inTreeFileName [2] inMissingListFileName [3] outTreeFileName");

	string inTreeFileName = argv[1];
	string inMissingListFileName =  argv[2];
	string outTreeFileName =  argv[3];	
	
	ofstream outTreesFile(outTreeFileName.c_str());
	vector<tree> allTrees = getStartingTreeVecFromFile(inTreeFileName);
	for (int t = 0; t < allTrees.size(); ++t) {
		removeTaxa(allTrees[t], inMissingListFileName);
		allTrees[t].output(outTreesFile);
	}
	outTreesFile.close();

	return 0;

}

void removeTaxa(tree& inTree,  const string& missingListFile){
	
	ifstream m_in(missingListFile.c_str());
	vector<string> names;
	putFileIntoVectorStringArray(m_in,names);

	for (int i=0; i<names.size();i++){
		tree::nodeP thisNode = inTree.findNodeByName(names[i]);
		if (thisNode) {
			if ((thisNode->father()->isRoot()) && (inTree.getRoot()->getNumberOfSons() == 2))
			{
				//in case the tree was rooted and the removed leaf was one of the root' sons:
				//we have to remove the root and reroot the tree at the second root son
				tree::nodeP pRoot = inTree.getRoot();
				tree::nodeP otherSonOfRoot;
				if (thisNode == pRoot->getSon(0))
					otherSonOfRoot = pRoot->getSon(1);
				else
					otherSonOfRoot = pRoot->getSon(0);

				tree newTree;
				newTree.createRootNode();
				newTree.getRoot()->setName(otherSonOfRoot->name());
				newTree.recursiveBuildTree(newTree.getRoot(),otherSonOfRoot->getSon(0));
				newTree.recursiveBuildTree(newTree.getRoot(),otherSonOfRoot->getSon(1));
				inTree = newTree;
			}
			else
                inTree.removeLeaf(thisNode);
		}
		else
		{
			cerr<<"Error in removeTaxa. Cannot find: "<<names[i]<<" in tree"<<endl;
			//errorMsg::reportError(err);
		}
	}
}

