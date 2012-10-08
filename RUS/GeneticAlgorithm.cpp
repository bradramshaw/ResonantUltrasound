#include "stdafx.h"
#include "GeneticAlgorithm.h"

static const int nVars = 9; // number of variables, F, dF, etc...
static const int nParams = 10;  //total parameters including chiSq

static const int nThreads = 4;
double totalTime = 0;  // for timing things
LARGE_INTEGER freq;
_CrtMemState s1,s2,s3;

GeneticAlgorithm::GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int order, double xHL, double yHL, double zHL, double density){	
	QueryPerformanceFrequency(&freq);
	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	

	_order = order;
	_R = 3 * (order+1) * (order+2) * (order+3) / 6;

	_xHL = xHL;
	_yHL = yHL;
	_zHL = zHL;

	_density = density;
	_basisPop = new int[8];
	for(int i = 0; i < 8; i++){
		_basisPop[i] = 0;
	}
	
	_basis = createBasis(order, _basisPop);

	initialiseMatrices();


	initializeParameters(dataSet, dataSetLength, nPopulation, scaleFactor, crossingProbability);	

	
	
}

GeneticAlgorithm::~GeneticAlgorithm(void){
	delete [] _populationParametersNew;
	delete [] _populationParametersOld;
	delete [] _basis;
}

void GeneticAlgorithm::initializeRandomNumberGenerators(){
	SYSTEMTIME t;
	GetLocalTime(&t);
	unsigned int max = _nPopulation - 1;

	vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
}

void GeneticAlgorithm::initialiseMatrices(){
	_emat = calcEmat(_R, _basis);
	_gradientCalcs = calcGradient(_R,_basis);
	return;
}

void GeneticAlgorithm::initializeParameters(double* dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability){

	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_dataSetLength = dataSetLength;
	_dataSet = dataSet;

	_residualArray = new double*[nThreads];
	_paramArray = new double*[nThreads];
	for(int i = 0; i < nThreads; i++){
	_residualArray[i] = new double[_dataSetLength];
	_paramArray[i] = new double[nParams];
	}
	
	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	for(int i  = 0; i < _nPopulation; i++){

		/*_populationParametersOld[i].c11 = (randomDouble(0,1))*pow(10,9);
		_populationParametersOld[i].c44 =  (randomDouble(0,1))*pow(10,9);*/

	/*	_populationParametersOld[i].c11 = (randomDouble(0,1))*pow(10,9);
		_populationParametersOld[i].c22 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c33 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c44 = (randomDouble(0,1))*pow(10,9);
		_populationParametersOld[i].c55 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c66 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c12 = (randomDouble(0,1))*pow(10,9);
		_populationParametersOld[i].c13 = _populationParametersOld[i].c12;
		_populationParametersOld[i].c23 = _populationParametersOld[i].c12;*/
		
		_populationParametersOld[i].c11 = 0.52296e+9;
		_populationParametersOld[i].c22 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c33 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c44 = 0.16288e+9;
		_populationParametersOld[i].c55 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c66 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c12 = 0.19721e+9;
		_populationParametersOld[i].c13 = _populationParametersOld[i].c12;
		_populationParametersOld[i].c23 = _populationParametersOld[i].c12;
	
		_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i],0);
	
	}

		_minimumParameters.c11 = 1;
		_minimumParameters.c22 = 1;
		_minimumParameters.c33 = 1;
		_minimumParameters.c44 = 1;
		_minimumParameters.c55 = 1;
		_minimumParameters.c66 = 1;
		_minimumParameters.c12 = 1;
		_minimumParameters.c13 = 1;
		_minimumParameters.c23 = 1;

		_minimumParameters.chiSq = INFINITE;
}

void GeneticAlgorithm::resetParameters(int nPopulation, double scaleFactor, double crossingProbability){
	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_nPopulation = nPopulation;
	delete [] _populationParametersOld;
	delete [] _populationParametersNew;

	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);

	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].c11 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c22 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c33 = _populationParametersOld[i].c11;
		_populationParametersOld[i].c44 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c55 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c66 = _populationParametersOld[i].c44;
		_populationParametersOld[i].c12 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c13 = _populationParametersOld[i].c12;
		_populationParametersOld[i].c23 = _populationParametersOld[i].c12;
		_populationParametersOld[i].chiSq = 1;
		_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i],0);
	}

	//delete _integerDistribution;
	delete [] ints1;
	delete [] ints2;
	delete [] ints3;

	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
	
}

double GeneticAlgorithm::randomDouble(double min, double max){				//change this to return full aray: faster
	double randNum;
	 vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &randNum, min, max);
	return randNum;
}

void GeneticAlgorithm::calculateMinimum(){
	for(int i  = 0; i < _nPopulation; i++){
		if(_populationParametersOld[i].chiSq < _minimumParameters.chiSq){
			_minimumParameters.c11 = _populationParametersOld[i].c11;
			_minimumParameters.c22 = _populationParametersOld[i].c22;
			_minimumParameters.c33 = _populationParametersOld[i].c33;
			_minimumParameters.c44 = _populationParametersOld[i].c44;
			_minimumParameters.c55 = _populationParametersOld[i].c55;
			_minimumParameters.c66 = _populationParametersOld[i].c66;
			_minimumParameters.c12 = _populationParametersOld[i].c12;
			_minimumParameters.c13 = _populationParametersOld[i].c13;
			_minimumParameters.c23 = _populationParametersOld[i].c23;			

			_minimumParameters.chiSq = _populationParametersOld[i].chiSq;
		}
	}
}

void GeneticAlgorithm::printMinimumParameters(){
	std::cout<<"c11: "<<_minimumParameters.c11<<" "<<"c22: "<<_minimumParameters.c22<<" "<<"c33: "<<_minimumParameters.c33<<std::endl<<"c44: "<<_minimumParameters.c44<<" "<<"c55: "<<_minimumParameters.c55<<" "<<"c66: "<<_minimumParameters.c66<<std::endl<<"c12: "<<_minimumParameters.c12<<" "<<"c13: "<<_minimumParameters.c13<<" "<<"c23: "<<_minimumParameters.c23<<std::endl<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
		std::ofstream out;
		out.open("output.dat");
		out.precision(15);
		out<<_minimumParameters.c11<<'\t'<<_minimumParameters.c22<<'\t'<<_minimumParameters.c33<<'\t'<<_minimumParameters.c44<<'\t'<<_minimumParameters.c55<<'\t'<<_minimumParameters.c66<<'\t'<<_minimumParameters.c12<<'\t'<<_minimumParameters.c13<<'\t'<<_minimumParameters.c23;
		out.close();
}

double GeneticAlgorithm::calculateResidual(Parameters::fitParameters * parameters, int threadID){
	
	double total = 0;	
	double * paramPointer = &(parameters->c11);	

	for(int i = 0; i < nParams; i++){
		_paramArray[threadID][i] = *paramPointer;
	
		paramPointer++;
	}
		
	double * frequencies = calculateFrequencies(_paramArray[threadID]);
	 

	
#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){

		_residualArray[threadID][i] = frequencies[i] - _dataSet[i];
		
	}	

#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){
		total += _residualArray[threadID][i]*_residualArray[threadID][i];
	}

	return total;
}

double * GeneticAlgorithm::calculateFrequencies(double * parameters){	
	double * frequencies = (double*) malloc(sizeof(double)*( _R - 6)) ;
	
	double **** ctens = initElasticConstants(parameters);


	LARGE_INTEGER time1,time2,freq;  // stores times and CPU frequency for profiling
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&time1);

	double * emat = new double[_R*_R];
	memcpy(emat, _emat, sizeof(double)*_R*_R);

	
	//QueryPerformanceCounter(&time2);
	//std::cout<<"Time to memcpy: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl;

	////	for(int i = 0; i < _R; i++){
	////	for(int j = 0; j < _R; j++){
	////		std::cout<<emat[i+_R*j]<<" ";
	////	}
	////	std::cout<<std::endl;
	////}
	////std::cout<<std::endl;
	//
	///*double * emat = calcEmat(_R, _basis);*/

	/*
	double * gmat = calcGmat(_R, _basis, ctens);*/
	

	//for(int i = 0; i < _R; i++){
	//	for(int j = 0; j < _R; j++){
	//		std::cout<<gmat[i+_R*j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//std::cout<<std::endl;

	double * gmat = calcGmat(_R, _basis, ctens, _gradientCalcs);
	/*		QueryPerformanceCounter(&time1);*/
	double * temp = calcEigs(_R, emat, gmat);
	//	QueryPerformanceCounter(&time2);
	//std::cout<<"Time to memcpy: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl;

	for(int i  = 6; i < _R; i++){
		frequencies[i-6] = (sqrt(temp[i]))/(2*3.1415926535897*pow(10,5));
	
	}
	

	delete [] emat;
	delete [] gmat;
	delete [] temp;
	
	return frequencies;
}

double **** GeneticAlgorithm::initElasticConstants(double * parameters){  // I shouldn't be storing this as a 4d array but it's easy this way
	double **** ctens;
		ctens = new double***[3];
		for(int i = 0; i <3; i++){
			ctens[i] = new double**[3];
				for(int j = 0; j<3; j++){
					ctens[i][j] = new double*[3];
					for(int k = 0; k<3; k++){
						ctens[i][j][k] = new double[3];
					}
				}
		}


		for(int i = 0; i < 3; i++){ //innitialize all entries to zero. this should be incorporated above obviously...
			for(int j = 0; j<3; j++){
				for(int k = 0; k<3; k++){
					for(int l = 0; l<3; l++){
						ctens[i][j][k][l] = 0;
					}
				}
			}
		}
		// note that the following takes care of symmetry, ie c_yzyz = c_zyzy = c_zyyz = c_yzzy
		ctens[0][0][0][0] = parameters[0];
	    ctens[0][0][1][1] = parameters[6];
		ctens[1][1][0][0] = parameters[6];
		ctens[0][0][2][2] = parameters[7];
		ctens[2][2][0][0] = parameters[7];
		ctens[1][1][1][1] = parameters[1];
		ctens[2][2][2][2] = parameters[2];
		ctens[2][2][1][1] = parameters[8];
		ctens[1][1][2][2] = parameters[8];

		ctens[1][2][1][2] = parameters[3];
		ctens[2][1][2][1] = parameters[3];
		ctens[2][1][1][2] = parameters[3];
		ctens[1][2][2][1] = parameters[3];
		
		ctens[2][0][2][0] = parameters[4];
		ctens[0][2][2][0] = parameters[4];
		ctens[0][2][0][2] = parameters[4];
		ctens[2][0][0][2] = parameters[4];

		ctens[0][1][0][1] = parameters[5];
		ctens[1][0][0][1] = parameters[5];
		ctens[0][1][1][0] = parameters[5];
		ctens[1][0][1][0] = parameters[5];

		return ctens;
}

double * GeneticAlgorithm::calcEigs(int R, double * emat, double * gmat){
		
		lapack_int ch0, ch, ch2;
		
	/*
		ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U',R, emat, R);*/
		int address = 0;
		for(int i = 0; i < 8; i++){
		
			ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', _basisPop[i], &emat[address], R);
			address += (_basisPop[i])*(R+1);
			
		}

		address = 0;
		for(int i = 0; i < 8; i++){
			ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U',  _basisPop[i], &gmat[address], R, &emat[address],R);
			address += (_basisPop[i])*(R+1);
		}
	/*	ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U', R, gmat, R, emat,R);*/

		double * w = (double*) malloc(sizeof(double)*R); // this will store eigenvalues

		address = 0;
		int position = 0;
		for(int i = 0; i < 8; i++){
			ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', _basisPop[i], &gmat[address], R, &w[position]);
			address += (_basisPop[i])*(R+1);
			position += _basisPop[i];
		}

		qsort(w,R,sizeof(double), dComp);
				
		//double * w = (double*) malloc(sizeof(double)*R); // this will store eigenvalues
		//ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', R, gmat, R, w); //computes eigenvalues and stores in "w". Eigenvectors are optional, but not computed here (Need to take advantage of block-diagonalization)
	
		return w;
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
	
	
	///*double * gmat;*/ // potential energy matrix. This is more complicated beacuse it depends on gradients 
		//double * gmat = new double[R*R]; // same size of course. 
	
		////again, this is symmetric, so only calculate the upper half part and just duplicate
		//for(int i = 0; i < R; i++){
		//	for(int j = i; j < R; j++){ 
		//		
		//		double tempSum = 0; // this is a sum of many terms; this is the storage variable
		//		for(int k = 0; k<3; k++){
		//			for(int l = 0; l < 3; l++){ // the elastic tensor is 4 dimensional. two coordinates come from the displacement directions the basis functions belong to, and the other two are the directions the derivatives are beign taken in. 
		//				tempSum += ctens[bFunctions[i].coord][k][bFunctions[j].coord][l]*integrateGradBasis(&bFunctions[i],&bFunctions[j],k,l, _xHL, _yHL, _zHL); //this just stupidly tries all of them, even though many are zero.
		//			}
		//				
		//		}
		//		gmat[i*R+j] = tempSum; // set upper and lower part of the matrix to the sum of the components. 
		//		gmat[j*R+i] = tempSum;
		//	}
		//}
	
		//return gmat;
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



void GeneticAlgorithm::calculateNewGenerations(int nGenerations){
	double averageTime = 0;

	for(int i = 0; i < nGenerations; i++){


	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints1, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints2, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints3, 0, _nPopulation );
	 				
		for(int j = 0; j < _nPopulation; j++){

		
			double * pointerToOldVariable = &_populationParametersOld[j].c11;
			double * pointerToNewVariable = &_populationParametersNew[j].c11;
			double * pointerTog1Variable = &_populationParametersOld[ints1[j]].c11;
			double * pointerTog2Variable = &_populationParametersOld[ints2[j]].c11;
			double * pointerTog3Variable = &_populationParametersOld[ints3[j]].c11;
		
	//could be optimzed for vector arithmetic
			for(int k = 0; k < nVars; k++){
				double p = randomDouble(0,1);
				if(p > _crossingProbability){
					*pointerToNewVariable = *pointerToOldVariable;
				}
				else{
					*pointerToNewVariable = *pointerTog1Variable + _scaleFactor*(*pointerTog2Variable - *pointerTog3Variable);
					if((k == 8) || (k == 9))
						if(*pointerToNewVariable > 1.0)
							*pointerToNewVariable = *pointerToNewVariable - 1.0;
						if(*pointerToNewVariable < -1.0)
							*pointerToNewVariable = *pointerToNewVariable + 1.0;
				}
				pointerToNewVariable++;
				pointerToOldVariable++;
				pointerTog1Variable++;
				pointerTog2Variable++;
				pointerTog3Variable++;
			}		
		}	

		totalTime = 0;

		HANDLE threadEvents[nThreads];
		Parameters::arrayBounds threadBounds[nThreads];
		threadContents threadContents[nThreads];
	
		for(int m = 0; m<nThreads; m++){		

			 threadEvents[m] = CreateEvent(NULL, FALSE, FALSE, NULL);
			int nPopulationPerThread =  _nPopulation/nThreads;
			if(m != (nThreads-1)){
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1;
			}
			else{
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1 + _nPopulation%nThreads;
			}		
			threadBounds[m].handle = threadEvents[m];
			threadBounds[m].time = 0;
			threadBounds[m].threadID = m;
			threadContents[m].arrayBounds = threadBounds[m];
			threadContents[m].pThis = this;

			AfxBeginThread(startResidualThread, (LPVOID) &threadContents[m]);		
		}
	
		//_CrtMemCheckpoint( &s2 );		
//		_CrtMemDumpStatistics( &s2 );
	//	std::cout<<s2.lTotalCount<<std::endl;

		WaitForMultipleObjects(nThreads,threadEvents,TRUE,INFINITE);	


		totalTime = threadContents[0].arrayBounds.time+threadContents[1].arrayBounds.time+threadContents[2].arrayBounds.time+threadContents[3].arrayBounds.time;
		averageTime +=totalTime;
	//	calculateMinimum();
	//	exportChiSq();
	}
	std::cout<<"Average time per generation per thread is: "<<averageTime/(4*nGenerations)<<"ms"<<std::endl<<std::endl;
	
}

UINT GeneticAlgorithm::startResidualThread(LPVOID param){
	threadContents * contents = (threadContents*) param;
	contents->pThis->residualCalculatingThread(&(contents->arrayBounds));
	return 0;
}

void GeneticAlgorithm::residualCalculatingThread(Parameters::arrayBounds * arrayBounds){
	LARGE_INTEGER time1,time2;
	#pragma ivdep
	for(int i = arrayBounds->start; i <= arrayBounds->end; i++){
		QueryPerformanceCounter(&time1);
		_populationParametersNew[i].chiSq = calculateResidual(&_populationParametersNew[i],arrayBounds->threadID);
		QueryPerformanceCounter(&time2);
		arrayBounds->time += 1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart);
			if(_populationParametersNew[i].chiSq < _populationParametersOld[i].chiSq){
				_populationParametersOld[i] = _populationParametersNew[i];
			}
	}
	SetEvent(arrayBounds->handle);
	return;
}

void GeneticAlgorithm::exportChiSq(){
    	std::ofstream out;
		out.open("chiSq.dat",std::ios_base::app);
		out.precision(15);
		out<<_minimumParameters.chiSq<<"\t";
		out.close();
}