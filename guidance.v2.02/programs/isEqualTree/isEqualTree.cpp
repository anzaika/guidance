#include <string>
#include <iostream>
#include "tree.h"
#include "treeUtil.h"

using namespace std;

int main(int argc, char** argv) {
	
	string treeFile1 = argv[1];
	string treeFile2 = argv[2];
	tree t1(treeFile1);
	tree t2(treeFile2);
	if (sameTreeTolopogy(t1,t2)) {
		cout << "1" << endl;
		return 1;
	}	else {
		cout << "0" << endl;
		return 0;
	}
	cout << "ERROR" << endl;
	return 2;
}
