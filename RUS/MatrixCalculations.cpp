#include "stdafx.h"
#include "GeneticAlgorithm.h"

// potential energy gmat, kintic emat, and calculate eigenvalues. Also contains a comparator for doulbes (for qsort). calcGradient calcualtes the part of the potential
//energy that only needs the dimensions (then gmat does the sum over the elastic tensor). The other two utility functions are integrateBasis and integrateGradBasis, used for the
//energy matrix calculations.

double * GeneticAlgorithm::calcEigs(int R, double * emat, double * gmat){
		
		lapack_int ch0, ch, ch2;	

		int address = 0;
		int addresses[8];
		int position = 0;
		int positions[8];
		for(int i = 0; i < 8; i ++){
			addresses[i] = address;
			positions[i] = position;
			address += (_basisPop[i])*(R+1);
			position += _basisPop[i];
		}
		
		
		/*#pragma ivdep
		for(int i = 0; i < 8; i++){
		
			ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', _basisPop[i], &emat[addresses[i]], R);

			
		}*/

	
	
		for(int i = 0; i < 8; i++){
	
			ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U',  _basisPop[i], &gmat[addresses[i]], R, &emat[addresses[i]],R);
		
	
		}

		double * w = (double*) malloc(sizeof(double)*R); // this will store eigenvalues



	
		for(int i = 0; i < 8; i++){
			ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', _basisPop[i], &gmat[addresses[i]], R, &w[positions[i]]);
		}

		qsort(w,R,sizeof(double), dComp);
	
		return w;
}


double * GeneticAlgorithm::calcGradient(int R, Basis::basisFunction * bFunctions){
	int elements = 0;
	for(int i = 0; i < 8; i++){
		elements += (_basisPop[i])*(_basisPop[i]);
	}

	double * gradientCalcs = new double[elements*9];

		int address = 0;
		int basisTotal = 0;
		int gradIndex = 0;
	
		//again, this is symmetric, so only calculate the upper half part and just duplicate
	for (int bN = 0; bN < 8; bN++){
		for(int i = 0; i < _basisPop[bN]; i++){
			for(int j = i; j < _basisPop[bN]; j++){ 
				
				 // this is a sum of many terms; this is the storage variable
				for(int k = 0; k<3; k++){
					for(int l = 0; l < 3; l++){ // the elastic tensor is 4 dimensional. two coordinates come from the displacement directions the basis functions belong to, and the other two are the directions the derivatives are beign taken in. 
						gradientCalcs[gradIndex]= integrateGradBasis(&bFunctions[basisTotal + i],&bFunctions[basisTotal + j],k,l, _xHL, _yHL, _zHL); //this just stupidly tries all of them, even though many are zero.
						gradIndex++;
					}
						
				}
			

			}
		}
		address += (_basisPop[bN])*(R+1);
		basisTotal += _basisPop[bN];
	}

	return gradientCalcs;
}

double * GeneticAlgorithm::calcGmat(int R, Basis::basisFunction * bFunctions, double **** ctens, double * gradientCalcs){
	
		/*double * gmat;*/ // potential energy matrix. This is more complicated beacuse it depends on gradients 
		double * gmat = new double[R*R]; // same size of course. 
		int address = 0;
		int basisTotal = 0;
		int gradIndex = 0;
		//again, this is symmetric, so only calculate the upper half part and just duplicate
	for (int bN = 0; bN < 8; bN++){
		for(int i = 0; i < _basisPop[bN]; i++){
			for(int j = i; j < _basisPop[bN]; j++){ 
				
				double tempSum = 0; // this is a sum of many terms; this is the storage variable
				for(int k = 0; k<3; k++){
					for(int l = 0; l < 3; l++){ // the elastic tensor is 4 dimensional. two coordinates come from the displacement directions the basis functions belong to, and the other two are the directions the derivatives are beign taken in. 
						//tempSum += ctens[bFunctions[basisTotal + i].coord][k][bFunctions[basisTotal + j].coord][l]*integrateGradBasis(&bFunctions[basisTotal + i],&bFunctions[basisTotal + j],k,l, _xHL, _yHL, _zHL); //this just stupidly tries all of them, even though many are zero.
						tempSum += ctens[bFunctions[basisTotal + i].coord][k][bFunctions[basisTotal + j].coord][l]*gradientCalcs[gradIndex]; //this just stupidly tries all of them, even though many are zero.

						gradIndex++;
					}
						
				}
				gmat[address + i*R + j] = tempSum; // set upper and lower part of the matrix to the sum of the components. 

		//		gmat[j*R+i] = tempSum;
			}
		}
		address += (_basisPop[bN])*(R+1);
		basisTotal += _basisPop[bN];
	}
	
		return gmat;		
	
}


double * GeneticAlgorithm::calcEmat(int R, Basis::basisFunction * bFunctions){

	double * emat; //pointer to the kinetic energy matrix. uses attrocious Fortran storage format (one-dimensional continuous array for a matrix...) for matrices because that's what Lapack wants.
	emat = new double[R*R]; // total size is of course the dimension squared.
	
		for(int i = 0; i < R; i++){  //calculates the kinetic energy matrix. 
			for(int j = i; j < R; j++){ // only calculates the upper-triangle, sine this is of course symmetric
				if(bFunctions[i].coord == bFunctions[j].coord){ // only basis functions belonging to the same coordinate contribute (each coordinate is a displacement direction, and each basis function part of an expansion of the displacement. kinetic energy obviously doesn't mix these). 
					
					double kinteticE = _density*integrateBasis(&bFunctions[i],&bFunctions[j], _xHL, _yHL, _zHL); // integrate these across the sample and multiply by the density to get the kinteic energy
					emat[R*i+j] =  kinteticE; // because this is a 1-D array, we have to be tricky in how we store it. i am storing both the upper and lower part
				//	emat[R*j+i] =  kinteticE; // lower part
				}
				else{
					emat[R*i+j] = 0; //should innitialize whole thing to zero and not deal with these cases explicitly. 
				//	emat[R*j+i] = 0;
				}
			}
		}
			
	
		return emat;
}


int GeneticAlgorithm::dComp(const void *a, const void *b){
	double diff = ( *(double*)a - *(double*)b );
	if(diff < 0)
		return -1;
	else return 1;
	
}


double GeneticAlgorithm::integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax){ // integrates two basis functions together. 
	double intVal = 8; // 2^3, because we integrate half.
	
	if( ( ( (*b1).xk + (*b2).xk)%2 == 1 )|| ( ( (*b1).yl + (*b2).yl)%2 == 1 ) || ( ( (*b1).zm + (*b2).zm)%2 == 1 ) ){ // if the sum of the powers for any coordinate are odd, then integral is zero
		intVal *= 0;
	}
	else{ // simple arithmetic for integration. obviously being a parallelapiped helps here.
		intVal *= ((1/( (double)(*b1).xk + (double)(*b2).xk + 1))*pow(xmax, (double)(*b1).xk + (double)(*b2).xk + 1))*((1/( (double)(*b1).yl + (double)(*b2).yl + 1))*pow(ymax, (double)(*b1).yl + (double)(*b2).yl + 1))*((1/( (double)(*b1).zm + (double)(*b2).zm + 1))*pow(zmax, (double)(*b1).zm + (double)(*b2).zm + 1));
	}
	return intVal;
}

double GeneticAlgorithm::integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax){ //integration for potential energy, where derivatives are taken
	double intVal = 8; //2^3, because integrate halves. 

	double powers[3] = {(double)(*b1).xk + (double)(*b2).xk,   (double)(*b1).yl + (double)(*b2).yl,   (double)(*b1).zm + (double)(*b2).zm}; // loads the total powers for each of x, y, and z
	
    int* ad1 = (int*) b1; // basis struct contains only ints, so we can use pointer arithemtic. this points to first integer in the struct (the x power)
	int* ad2 = (int*) b2; //same for second basis function

	ad1 += d1; //move pointers to whichever coordinate the derviative is being taken with respect to. 
	ad2 += d2;

	intVal *= ((double) *ad1)*((double) *ad2); //multiply by the coefficients we get from taking a derivative.

	powers[d1]--; // decrement the powers on the coordinates that had the derivatives taken.
	powers[d2]--;

	if((powers[0] <0 )||(powers[1] <0 )||(powers[2] <0 )){ // if we took derivatives of a coordinate that was 0th order (eg x^0 = 1) then we get zero for the integral
		return 0;
	}

	if((( (int) powers[0]%2 == 1) ||( (int) powers[1]%2 == 1)||( (int) powers[2]%2 == 1)) ){ //also get zero for integrating odd powers
		return 0;
	}

	intVal *= (1/(powers[0]+1)*pow(xmax,powers[0]+1))*(1/(powers[1]+1)*pow(ymax,powers[1]+1))*(1/(powers[2]+1)*pow(zmax,powers[2]+1)); //otherwise simple arithmetic
	
	return intVal;
}

