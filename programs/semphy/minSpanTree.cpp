// $Id: minSpanTree.cpp 409 2005-06-28 13:12:24Z ninio $

#include "minSpanTree.h"
#include "logFile.h"
#include <queue>
using namespace std;

rearrangeTree::pairSet minSpanTree::span_tree(const VVdouble &weight,
												double* out_score){

	MDOUBLE tmpWeight;
	int numOfNodes = weight.size();
	priority_queue<edgePair> edges;
	vector<bool> notDone(numOfNodes,true);
	int nodesLeft = numOfNodes;
	rearrangeTree::pairSet edgeList;
	double score =0.0;
	int new_node_id=0;

	//  cout << weight<<endl;
	notDone[new_node_id]=false;
	--nodesLeft;
	for (int i=0;i<numOfNodes;i++) {
		if (notDone[i] && i != new_node_id) {
			tmpWeight = (new_node_id<i?weight[new_node_id][i]:weight[i][new_node_id]);
			edges.push(edgePair(tmpWeight, rearrangeTree::intPair(new_node_id,i)));
		}
	}

	while (nodesLeft>0) {
		edgePair best = edges.top(); 
		LOG(50,<<"trying edge ["<<best.second.first<<","<<best.second.second<<"]="<<best.first<<endl);
		edges.pop();		// remove top
		if (!notDone[best.second.first]&&!notDone[best.second.second]) //not one of the sides is new
		continue;
		else {
		  LOG(50,<<"accepted - adding to score"<<endl);
			score += best.first;
			//	cout << best.first<<"\t"<< best.second.first<<"\t"<<best.second.second<<endl;
			edgeList.insert(best.second);
			if (notDone[best.second.first]) {
				new_node_id=best.second.first;
			} 
			else {
				new_node_id=best.second.second;
			}
			notDone[new_node_id]=false;
			--nodesLeft;
			for (int i=0;i<numOfNodes;i++) {
				if (notDone[i]){
					tmpWeight = (new_node_id<i?weight[new_node_id][i]:weight[i][new_node_id]);
					edges.push(edgePair(tmpWeight,rearrangeTree::intPair(new_node_id,i)));
				}
			}
		}
	}
	*out_score=score;
	return edgeList;
}


