// $Id: weightVectors.cpp 409 2005-06-28 13:12:24Z ninio $

#include "weightVectors.h"
#include "someUtil.h"
#include "errorMsg.h"
#include "getRandomWeights.h"
#include <cassert>

// this is the format of the weight vectors

//5 2
//
//3   56 4 31 6
//653 4  4 0  0.4
//
//In the first row, the first number is the number of positions. This should
//be compatible with the sequence lengths.
//The second number is the nubmer of repeats (samples).
//Note that the weights can be floats or doubles.


weightVectors::weightVectors(istream& infile) {
	read(infile);
}

weightVectors::weightVectors(const int length, const int size) {
	checkLengthAndSize(length,size);
	_weights.resize(size);
	for (int i=0; i < _weights.size();++i) {
		_weights[i].resize(length,1.0);
	}
}

void weightVectors::randomizeBPweights(){
	for (int i=0; i < _weights.size();++i) {
		getRandomWeights::standardBPWeights(_weights[i]);
	}
}

void weightVectors::checkLengthAndSize(const int length, const int size) {
if (length <=0) {
		errorMsg::reportError("problem reading weight vectors. The length should be a positive integer");
	}
	if (size <=0) {
		errorMsg::reportError("problem reading weight vectors. The number of weight vectors should be a positive integer");
	}

}
int weightVectors::read(istream &infile) {
	// reading the first line
	int length; infile >> length;
	int size; infile >> size;
	checkLengthAndSize(length,size);
	_weights.resize(size);
	// reading the vectors
	int i;
	for (i=0; i < size; ++i) {
		_weights[i].resize(length);
		for (int j=0; j < length; ++j) {
			if (infile.eof()) {
				string errMsg = "problem reading weight vectors.\n";
				errMsg+="expected number of vectors was ";
				errMsg+=int2string(size)+".\n";
				errMsg+="expected length of each vector was ";
				errMsg+=int2string(length)+"\n";
				errMsg+="Read "+int2string(i)+" full vectors and "+int2string(j)+" numbers\n";
				errorMsg::reportError(errMsg);
			}
			infile>>_weights[i][j];
		}
	}
	return i;
}

