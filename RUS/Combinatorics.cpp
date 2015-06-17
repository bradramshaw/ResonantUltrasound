#include "stdafx.h"
#include "GeneticAlgorithm.h"

void GeneticAlgorithm::shiftInd(long * indV,long shiftIndex, long nMiss, long nFreq, bool * maxFlag){
 
	if((shiftIndex == 0)&& (indV[shiftIndex] == nFreq - (nMiss - shiftIndex))){
		return;
	}	

	if(indV[shiftIndex] == nFreq -(nMiss - shiftIndex)){
		*maxFlag = true;
		shiftInd(indV, shiftIndex - 1, nMiss, nFreq, maxFlag);
	}
	
	else{
	*maxFlag = false;
	
	indV[shiftIndex]++;

	for(int i = shiftIndex +1; i < nMiss; i ++){
		indV[i] = indV[shiftIndex]+(i-shiftIndex );
	}

	return;
	}
}



long GeneticAlgorithm::nCombs(int N, int k){

	long result = (long) N;
	long kfac = (long) k;

	if(N == 0){
		return 0;
	}
	if(k == 0){
		return 1;
	}

	for(int i = N-1; i > (N-k); i--){
		result *= i;
	}
	for(int i = k-1; i >0; i--){
		kfac *= i;
	}

	result /= kfac;

	return result;
}