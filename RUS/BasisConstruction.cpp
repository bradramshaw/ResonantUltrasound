#include "stdafx.h"
#include "GeneticAlgorithm.h"

// create the basis, cacluate the parity, of each function, and a comparator for parity.

Basis::basisFunction * GeneticAlgorithm::createBasis(int order, int* basisPop){
	int R = 3 * (order+1) * (order+2) * (order+3) / 6;
	Basis::basisFunction * bFunctions = (Basis::basisFunction *) malloc(R * sizeof(Basis::basisFunction)); // allocates memory for the basis functions, of which there are R
		int basisPoint = 0; //basis functions, say p_i, are repeated for the x, y, and z coordinates. So we have p_i for x, for y, and for z. basisPoint is the "i" index in this notation.

		for(int k = 0; k <= order ; k++){  // k, l, and m are the powers for x, y, and z. 
			for(int l = 0; l <= order; l++){ // maximum size for any of k, l, and m is the order (x^4 is highest x power for 4th order, for example)
				for(int m = 0; m <= order; m++){
					if(k + l + m <= order){ // the sum of k, l, and m must be less than the order ( x*y*z^2 is 4th order)
						for(int i = 0; i<3; i++){ // i is for each of the three coordinates that each basis function is repeated for. The function is the same for each.
							bFunctions[basisPoint+i].xk = k;
							bFunctions[basisPoint+i].yl = l;
							bFunctions[basisPoint+i].zm = m;

							bFunctions[basisPoint+i].coord = i; // whether this particular basis function belongs to x, y, or z.
							bFunctions[basisPoint+i].pnumb = parity(k,l,m, i, basisPop); //calls the partity function to determine the parity. see albert's paper.
						}
					basisPoint += 3; // skip ahead three in the array because basis functions come in triplets (one for each coordinate). 
					}
				}
			}
		}
		qsort(bFunctions, R, sizeof(Basis::basisFunction), GeneticAlgorithm::compPnumb);
		return bFunctions;
}

int GeneticAlgorithm::parity(int k, int l, int m, int i, int * basisPop){ // calcualtes the parity of the function (functions of different parity integrate to zero, so this block-diagonalizes things). 
	double pvec[3] = {0,0,0};

	if(i == 0){ // i is the displacement coordinate of the function, x y or z
		pvec[0] =  pow(-1,k+1) ; // see albert's "potato" paper for definitions. 
		pvec[1] =  pow(-1,l) ;
		pvec[2] =  pow(-1,m) ;
	}
	else if(i == 1){
		pvec[0] =  pow(-1,k) ;
		pvec[1] =  pow(-1,l+1) ;
		pvec[2] =  pow(-1,m) ;
	}
	else if(i == 2){
		pvec[0] =  pow(-1,k) ;
		pvec[1] =  pow(-1,l) ;
		pvec[2] =  pow(-1,m+1) ;
	}


	if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == 1)){ //this categorizes the parities in x,y and z into eight groups. 
		basisPop[0]++;
	return 1;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
		basisPop[1]++;
	return 2;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
		basisPop[2]++;
	return 3;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
		basisPop[3]++;
	return 4;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == 1)){
		basisPop[4]++;
	return 5;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
		basisPop[5]++;
	return 6;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
		basisPop[6]++;
	return 7;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
		basisPop[7]++;
	return 8;
	}
	else return -1;
}


int GeneticAlgorithm::compPnumb(const void* b1, const void* b2){ //compartor for basis function structs. just looks at the parity lable, from 1 to 8.
	Basis::basisFunction * bas1 = (Basis::basisFunction*) b1;
	Basis::basisFunction * bas2 = (Basis::basisFunction*) b2;

	return 	(*bas1).pnumb - (*bas2).pnumb;

}