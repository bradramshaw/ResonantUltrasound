// RUS.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
using namespace std;

// All of the non-zero elastic constants, in pascals. The notation is 1 xx, 2 yy, 3 zz, 4 yz, 5 xz, 6 xy. there are symmetries, and so for example c12 = c_xxyy = c_yyxx = c21
static const double c11 = 0.52296e+9; 
static const double c22 = 0.52296e+9;
static const double c33 = 0.52296e+9;
static const double c44 = 0.16288e+9;
static const double c55 = 0.16288e+9;
static const double c66 = 0.16288e+9;
static const double c12 = 0.19721e+9;
static const double c13 = 0.19721e+9;
static const double c23 = 0.19721e+9;



static const double density = 15467; // density of the material in grams/meter^3. All units are SI

static const double xHL = 0.0023724; // HALF length in the x direction, in meters.
static const double yHL = 0.001885; // HALF length in the y direction
static const double zHL = 0.00242095; // HALF length in the z direction.


//int parity(int k, int l, int m, int coord); // parity function looks at the symmetry of a basis function, ie x^2 * y * z^3. more below in the full function definition.
//int compPnumb(const void * b1, const void * b2); //comparitor for parity numbers. 
//double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax); // This integrates  a pair of basis functions within the limits specified. Note that this assumes a parallelapiped, and takes "half" dimensions as inputs
//double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax); // integrates basis functions after differentiating with respect to one coordinate in each basis function of the pair. 
//double **** initElasticConstants(); //initialize the full tensor. 
//Basis::basisFunction * createBasis(int order);
//double * calcEmat(int order,Basis::basisFunction * bFunctions);  // kinetic energy
//double * calcGmat(int order, Basis::basisFunction * bFunctions, double **** ctens); // elastic energy
//double * calcEigs(int order, double * emat, double * gmat);// eigenvalues (resonant frequencies squared, or maybe their inverse. anyway it's obvious)

int _tmain(int argc, _TCHAR* argv[]) //main function
{
	LARGE_INTEGER time1,time2,freq;  // stores times and CPU frequency for profiling
	QueryPerformanceFrequency(&freq);
	//double **** ctens = initElasticConstants(); // 4 dimensional elastic constant array. can probably be simpler (obviously)

	
		VSLStreamStatePtr stream;
		SYSTEMTIME t;
	    GetLocalTime(&t);
	    vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
	
		DataExtractor extractor("E:/Users/Brad/Documents/GitHub/ResonantUltrasound/RUS/IsotropicCuboid.dat");
		double * data = extractor.getDataArray();
		int nPoints = extractor.getNumberOfLines();
	
		int order; // will store the max order of the polynomials to use
		double scale, cross;
		cout << "Highest polynomial order? ";
		cin >> order;
		cout  << endl;
		int R = 3 * (order+1) * (order+2) * (order+3) / 6; // total dimension of the matrices 
		cout << "R = " << R<<endl; // output that dimension to the user
		cout << "Scale factor? ";
		cin >> scale;
		cout  << endl;

		cout << "Crossing Probability? ";
		cin >> cross;
		cout << endl;

		GeneticAlgorithm geneticAlgorithm(data, nPoints,  40, scale, cross, order, xHL, yHL, zHL, density);
	
		geneticAlgorithm.calculateMinimum();
		geneticAlgorithm.printMinimumParameters();	

	while(true){ // bad programming...

		int nGens;
		cout<<"Number of generations: ";
		cin>>nGens;
		cout<<endl;
				
		geneticAlgorithm.calculateNewGenerations(nGens);		
	

		geneticAlgorithm.calculateMinimum();
		
		geneticAlgorithm.printMinimumParameters();	

		//Basis::basisFunction * bFunctions = createBasis(order);
		//	 // there are eight possible parities, and the basis functions are labled as such. This sorts them using the comparator compPnumb				
		//double * emat = calcEmat(order, bFunctions);
	
		//double * gmat =  calcGmat(order,  bFunctions,  ctens);

		////QueryPerformanceCounter(&time1);
		////lapack_int ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', R, emat, R); // this performs a cholesky decomposition on the kinetic matrix. result is stored in emat (cholesky decomposition on A gives A = L L*, where L is lower triangular. Only Hermitian matrices need apply)
		////QueryPerformanceCounter(&time2);
		////cout<<"Time to cholesky: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<endl<<endl;

		////QueryPerformanceCounter(&time1);
		////lapack_int ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U', R, gmat, R, emat,R); // reduces a generalized symmetric-definite generalized eignevalue problem to a standard eigenvalue problem
		////// problem is currently of type A*z = λ*B*z, this changes it to type C*z = λ*z. Essentially gives Binv*A*z = λ*z, but more clever (or stable. see intel MKL library documentation). 
		////QueryPerformanceCounter(&time2);
		////cout<<"Time to invert and reduce: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<endl<<endl;

		////QueryPerformanceCounter(&time1);
		////double * w = (double*) malloc(sizeof(double)*R); // this will store eigenvalues
		////lapack_int ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', R, gmat, R, w); //computes eigenvalues and stores in "w". Eigenvectors are optional, but not computed here (Need to take advantage of block-diagonalization)
		////QueryPerformanceCounter(&time2);
		////cout<<"Time to solve eigenproblem: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<endl<<endl;	

		//double * eigs = calcEigs(order,  emat,  gmat);

		//std::ofstream out4;
		//out4.open("eigs.dat",std::ios_base::beg); //writes eigenvalues so that they can be opened in mathematica or something.
		//out4.precision(15);
	
		//for(int i = 0; i < R; i++){
		// out4<<eigs[i]<<"\t"; //outputs the eigenvalues to the file
		//}
		//out4.close();

		/*std::ofstream out3;
		out3.open("imat.dat",std::ios_base::beg);
		out3.precision(15);
		int index3 = 0;
		for(int i = 0; i < R; i++){
			for(int j = 0; j<R; j++){
			out3<<gmat[index3]<<"\t";
			index3++;
			}
			out3<<"\n";
		}
		out3.close();*/

		//delete[] emat; //cleanup
		//delete[] gmat;
		//delete[] eigs;
		//delete[] bFunctions;

		}
	return 0;
}


//double * calcEigs(int order, double * emat, double * gmat){
//
//		int R = 3 * (order+1) * (order+2) * (order+3) / 6;
//		lapack_int ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', R, emat, R);
//		lapack_int ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U', R, gmat, R, emat,R);
//		double * w = (double*) malloc(sizeof(double)*R); // this will store eigenvalues
//		lapack_int ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', R, gmat, R, w); //computes eigenvalues and stores in "w". Eigenvectors are optional, but not computed here (Need to take advantage of block-diagonalization)
//
//		return w;
//}

//
//double * calcGmat(int order, Basis::basisFunction * bFunctions, double **** ctens){
//		int R = 3 * (order+1) * (order+2) * (order+3) / 6;
//		double * gmat; // potential energy matrix. This is more complicated beacuse it depends on gradients 
//		gmat = new double[R*R]; // same size of course. 
//	
//		//again, this is symmetric, so only calculate the upper half part and just duplicate
//		for(int i = 0; i < R; i++){
//			for(int j = i; j < R; j++){ 
//				
//				double tempSum = 0; // this is a sum of many terms; this is the storage variable
//				for(int k = 0; k<3; k++){
//					for(int l = 0; l < 3; l++){ // the elastic tensor is 4 dimensional. two coordinates come from the displacement directions the basis functions belong to, and the other two are the directions the derivatives are beign taken in. 
//						tempSum += ctens[bFunctions[i].coord][k][bFunctions[j].coord][l]*integrateGradBasis(&bFunctions[i],&bFunctions[j],k,l, xHL, yHL, zHL); //this just stupidly tries all of them, even though many are zero.
//					}
//				}
//				gmat[i*R+j] = tempSum; // set upper and lower part of the matrix to the sum of the components. 
//				gmat[j*R+i] = tempSum;
//			}
//		}
//		return gmat;
//}



//double * calcEmat(int order, Basis::basisFunction * bFunctions){
//
//	int R = 3 * (order+1) * (order+2) * (order+3) / 6;
//	double * emat; //pointer to the kinetic energy matrix. uses attrocious Fortran storage format (one-dimensional continuous array for a matrix...) for matrices because that's what Lapack wants.
//		emat = new double[R*R]; // total size is of course the dimension squared.
//	
//		for(int i = 0; i < R; i++){  //calculates the kinetic energy matrix. 
//			for(int j = i; j < R; j++){ // only calculates the upper-triangle, sine this is of course symmetric
//				if(bFunctions[i].coord == bFunctions[j].coord){ // only basis functions belonging to the same coordinate contribute (each coordinate is a displacement direction, and each basis function part of an expansion of the displacement. kinetic energy obviously doesn't mix these). 
//					double kinteticE = density*integrateBasis(&bFunctions[i],&bFunctions[j], xHL, yHL, zHL); // integrate these across the sample and multiply by the density to get the kinteic energy
//					emat[R*i+j] =  kinteticE; // because this is a 1-D array, we have to be tricky in how we store it. i am storing both the upper and lower part
//					emat[R*j+i] =  kinteticE; // lower part
//				}
//				else{
//					emat[R*i+j] = 0; //should innitialize whole thing to zero and not deal with these cases explicitly. 
//					emat[R*j+i] = 0;
//				}
//			}
//		}
//		return emat;
//}


//Basis::basisFunction * createBasis(int order){
//	int R = 3 * (order+1) * (order+2) * (order+3) / 6;
//	Basis::basisFunction * bFunctions = (Basis::basisFunction *) malloc(R * sizeof(Basis::basisFunction)); // allocates memory for the basis functions, of which there are R
//		int basisPoint = 0; //basis functions, say p_i, are repeated for the x, y, and z coordinates. So we have p_i for x, for y, and for z. basisPoint is the "i" index in this notation.
//
//		for(int k = 0; k <= order ; k++){  // k, l, and m are the powers for x, y, and z. 
//			for(int l = 0; l <= order; l++){ // maximum size for any of k, l, and m is the order (x^4 is highest x power for 4th order, for example)
//				for(int m = 0; m <= order; m++){
//					if(k + l + m <= order){ // the sum of k, l, and m must be less than the order ( x*y*z^2 is 4th order)
//						for(int i = 0; i<3; i++){ // i is for each of the three coordinates that each basis function is repeated for. The function is the same for each.
//							bFunctions[basisPoint+i].xk = k;
//							bFunctions[basisPoint+i].yl = l;
//							bFunctions[basisPoint+i].zm = m;
//
//							bFunctions[basisPoint+i].coord = i; // whether this particular basis function belongs to x, y, or z.
//							bFunctions[basisPoint+i].pnumb = parity(k,l,m, i); //calls the partity function to determine the parity. see albert's paper.
//						}
//					basisPoint += 3; // skip ahead three in the array because basis functions come in triplets (one for each coordinate). 
//					}
//				}
//			}
//		}
//		qsort(bFunctions, R, sizeof(Basis::basisFunction), compPnumb);
//		return bFunctions;
//}


//int parity(int k, int l, int m, int i){ // calcualtes the parity of the function (functions of different parity integrate to zero, so this block-diagonalizes things). 
//	double pvec[3] = {0,0,0};
//
//	if(i == 0){ // i is the displacement coordinate of the function, x y or z
//		pvec[0] =  pow(-1,k+1) ; // see albert's "potato" paper for definitions. 
//		pvec[1] =  pow(-1,l) ;
//		pvec[2] =  pow(-1,m) ;
//	}
//	else if(i == 1){
//		pvec[0] =  pow(-1,k) ;
//		pvec[1] =  pow(-1,l+1) ;
//		pvec[2] =  pow(-1,m) ;
//	}
//	else if(i == 2){
//		pvec[0] =  pow(-1,k) ;
//		pvec[1] =  pow(-1,l) ;
//		pvec[2] =  pow(-1,m+1) ;
//	}
//
//
//	if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == 1)){ //this categorizes the parities in x,y and z into eight groups. 
//	return 1;
//	}
//	else if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
//	return 2;
//	}
//	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
//	return 3;
//	}
//	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
//	return 4;
//	}
//	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == 1)){
//	return 5;
//	}
//	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
//	return 6;
//	}
//	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
//	return 7;
//	}
//	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
//	return 8;
//	}
//	else return -1;
//}

//int compPnumb(const void* b1, const void* b2){ //compartor for basis function structs. just looks at the parity lable, from 1 to 8.
//	Basis::basisFunction * bas1 = (Basis::basisFunction*) b1;
//	Basis::basisFunction * bas2 = (Basis::basisFunction*) b2;
//
//	return 	(*bas1).pnumb - (*bas2).pnumb;
//
//}
//
//double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax){ // integrates two basis functions together. 
//	double intVal = 8; // 2^3, because we integrate half.
//
//	if( ( ( (*b1).xk + (*b2).xk)%2 == 1 )|| ( ( (*b1).yl + (*b2).yl)%2 == 1 ) || ( ( (*b1).zm + (*b2).zm)%2 == 1 ) ){ // if the sum of the powers for any coordinate are odd, then integral is zero
//		intVal *= 0;
//	}
//	else{ // simple arithmetic for integration. obviously being a parallelapiped helps here.
//		intVal *= ((1/( (double)(*b1).xk + (double)(*b2).xk + 1))*pow(xmax, (double)(*b1).xk + (double)(*b2).xk + 1))*((1/( (double)(*b1).yl + (double)(*b2).yl + 1))*pow(ymax, (double)(*b1).yl + (double)(*b2).yl + 1))*((1/( (double)(*b1).zm + (double)(*b2).zm + 1))*pow(zmax, (double)(*b1).zm + (double)(*b2).zm + 1));
//	}
//	return intVal;
//}
//
//double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax){ //integration for potential energy, where derivatives are taken
//	double intVal = 8; //2^3, because integrate halves. 
//
//	double powers[3] = {(double)(*b1).xk + (double)(*b2).xk,   (double)(*b1).yl + (double)(*b2).yl,   (double)(*b1).zm + (double)(*b2).zm}; // loads the total powers for each of x, y, and z
//	
//    int* ad1 = (int*) b1; // basis struct contains only ints, so we can use pointer arithemtic. this points to first integer in the struct (the x power)
//	int* ad2 = (int*) b2; //same for second basis function
//
//	ad1 += d1; //move pointers to whichever coordinate the derviative is being taken with respect to. 
//	ad2 += d2;
//
//	intVal *= ((double) *ad1)*((double) *ad2); //multiply by the coefficients we get from taking a derivative.
//
//	powers[d1]--; // decrement the powers on the coordinates that had the derivatives taken.
//	powers[d2]--;
//
//	if((powers[0] <0 )||(powers[1] <0 )||(powers[2] <0 )){ // if we took derivatives of a coordinate that was 0th order (eg x^0 = 1) then we get zero for the integral
//		return 0;
//	}
//
//	if((( (int) powers[0]%2 == 1) ||( (int) powers[1]%2 == 1)||( (int) powers[2]%2 == 1)) ){ //also get zero for integrating odd powers
//		return 0;
//	}
//
//	intVal *= (1/(powers[0]+1)*pow(xmax,powers[0]+1))*(1/(powers[1]+1)*pow(ymax,powers[1]+1))*(1/(powers[2]+1)*pow(zmax,powers[2]+1)); //otherwise simple arithmetic
//	
//	return intVal;
//}
//
//double **** initElasticConstants(){  // I shouldn't be storing this as a 4d array but it's easy this way
//	double **** ctens;
//		ctens = new double***[3];
//		for(int i = 0; i <3; i++){
//			ctens[i] = new double**[3];
//				for(int j = 0; j<3; j++){
//					ctens[i][j] = new double*[3];
//					for(int k = 0; k<3; k++){
//						ctens[i][j][k] = new double[3];
//					}
//				}
//		}
//
//
//		for(int i = 0; i < 3; i++){ //innitialize all entries to zero. this should be incorporated above obviously...
//			for(int j = 0; j<3; j++){
//				for(int k = 0; k<3; k++){
//					for(int l = 0; l<3; l++){
//						ctens[i][j][k][l] = 0;
//					}
//				}
//			}
//		}
//		// note that the following takes care of symmetry, ie c_yzyz = c_zyzy = c_zyyz = c_yzzy
//		ctens[0][0][0][0] = c11;
//	    ctens[0][0][1][1] = c12;
//		ctens[1][1][0][0] = c12;
//		ctens[0][0][2][2] = c13;
//		ctens[2][2][0][0] = c13;
//		ctens[1][1][1][1] = c22;
//		ctens[2][2][2][2] = c33;
//		ctens[2][2][1][1] = c23;
//		ctens[1][1][2][2] = c23;
//
//		ctens[1][2][1][2] = c44;
//		ctens[2][1][2][1] = c44;
//		ctens[2][1][1][2] = c44;
//		ctens[1][2][2][1] = c44;
//		
//		ctens[2][0][2][0] = c55;
//		ctens[0][2][2][0] = c55;
//		ctens[0][2][0][2] = c55;
//		ctens[2][0][0][2] = c55;
//
//		ctens[0][1][0][1] = c66;
//		ctens[1][0][0][1] = c66;
//		ctens[0][1][1][0] = c66;
//		ctens[1][0][1][0] = c66;
//
//		return ctens;
//}