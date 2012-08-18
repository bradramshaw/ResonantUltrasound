// RUS.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
using namespace std;

static const double c11 = 231e+9;
static const double c22 = 268e+9;
static const double c33 = 186e+9;
static const double c44 = 49e+9;
static const double c55 = 37e+9;
static const double c66 = 95e+9;
static const double c12 = 132e+9;
static const double c13 = 71e+9;
static const double c23 = 95e+9;

static const double density = 6394;

static const double xHL = 0.000515;
static const double yHL = 0.0006;
static const double zHL = 0.0001025;


int parity(int k, int l, int m, int coord);
int compPnumb(const void * b1, const void * b2);
double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax);
double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax);

int _tmain(int argc, _TCHAR* argv[])
{
	char b;
	while(true){
	
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


		for(int i = 0; i < 3; i++){
			for(int j = 0; j<3; j++){
				for(int k = 0; k<3; k++){
					for(int l = 0; l<3; l++){
						ctens[i][j][k][l] = 0;
					}
				}
			}
		}

		ctens[0][0][0][0] = c11;
	    ctens[0][0][1][1] = c12;
		ctens[1][1][0][0] = c12;
		ctens[0][0][2][2] = c13;
		ctens[2][2][0][0] = c13;
		ctens[1][1][1][1] = c22;
		ctens[2][2][2][2] = c33;
		ctens[2][2][1][1] = c23;
		ctens[1][1][2][2] = c23;

		ctens[1][2][1][2] = c44;
		ctens[2][1][2][1] = c44;
		ctens[2][1][1][2] = c44;
		ctens[1][2][2][1] = c44;
		
		ctens[2][0][2][0] = c55;
		ctens[0][2][2][0] = c55;
		ctens[0][2][0][2] = c55;
		ctens[2][0][0][2] = c55;

		ctens[0][1][0][1] = c66;
		ctens[1][0][0][1] = c66;
		ctens[0][1][1][0] = c66;
		ctens[1][0][1][0] = c66;

		
		int order;
		cout << "Highest polynomial order? ";
		cin >> order;
		cout  << endl;
		int R = 3 * (order+1) * (order+2) * (order+3) / 6;
		cout << "R = " << R<<endl;

		Basis::basisFunction * bFunctions = (Basis::basisFunction *) malloc(R * sizeof(Basis::basisFunction));
		int basisPoint = 0;

		for(int k = 0; k <= order ; k++){
			for(int l = 0; l <= order; l++){
				for(int m = 0; m <= order; m++){
					if(k + l + m <= order){
						for(int i = 0; i<3; i++){
							bFunctions[basisPoint+i].xk = k;
							bFunctions[basisPoint+i].yl = l;
							bFunctions[basisPoint+i].zm = m;

							bFunctions[basisPoint+i].coord = i;
							bFunctions[basisPoint+i].pnumb = parity(k,l,m, i);
						}
					basisPoint += 3;
					}
				}
			}
		}
		 
		qsort(bFunctions, R, sizeof(Basis::basisFunction), compPnumb);
				
	/*	double ** emat;
		emat = new double*[R];
		for(int i = 0; i < R; i++){
		emat[i] = new double[R];
		}
		
		for(int i = 0; i < R; i++){
			for(int j = 0; j < R; j++){
				if(bFunctions[i].coord == bFunctions[j].coord){
					emat[i][j] =  15467*integrateBasis(&bFunctions[i],&bFunctions[j], 0.0023724, 0.00242095, 0.001885);
				}
				else{
					emat[i][j] = 0;
				}
			}
		}*/
		
		double * emat;
		emat = new double[R*R];
	
		for(int i = 0; i < R; i++){
			for(int j = i; j < R; j++){
				if(bFunctions[i].coord == bFunctions[j].coord){
					double kinteticE = density*integrateBasis(&bFunctions[i],&bFunctions[j], xHL, yHL, zHL);
					emat[R*i+j] =  kinteticE;
					emat[R*j+i] =  kinteticE;
				}
				else{
					emat[R*i+j] = 0;
					emat[R*j+i] = 0;
				}
			}
		}
	/*	std::ofstream out;
		out.open("emat.dat",std::ios_base::beg);
		out.precision(15);
		int index = 0;
		for(int i = 0; i < R; i++){
			for(int j = 0; j<R; j++){
			out<<emat[index]<<"\t";
			index++;
			}
			out<<"\n";
		}
		out.close();*/
		
		double * gmat;
		gmat = new double[R*R];
	
		for(int i = 0; i < R; i++){
			for(int j = i; j < R; j++){
				
				double tempSum = 0;
				for(int k = 0; k<3; k++){
					for(int l = 0; l < 3; l++){
						tempSum += ctens[bFunctions[i].coord][k][bFunctions[j].coord][l]*integrateGradBasis(&bFunctions[i],&bFunctions[j],k,l, xHL, yHL, zHL);
					}
				}
				gmat[i*R+j] = tempSum;
				gmat[j*R+i] = tempSum;
			}
		}


		//std::ofstream out2;
		//out2.open("gmat.dat",std::ios_base::beg);
		//out2.precision(15);
		//int index2 = 0;
		//for(int i = 0; i < R; i++){
		//	for(int j = 0; j<R; j++){
		//	out2<<gmat[index2]<<"\t";
		//	index2++;
		//	}
		//	out2<<"\n";
		//}
		//out2.close();

		lapack_int ch0 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', R, emat, R);
		lapack_int ch = LAPACKE_dsygst(LAPACK_ROW_MAJOR, 1,'U', R, gmat, R, emat,R);
		double * w = (double*) malloc(sizeof(double)*R);
		lapack_int ch2 = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', R, gmat, R, w);

		std::ofstream out4;
		out4.open("eigs.dat",std::ios_base::beg);
		out4.precision(15);
	
		for(int i = 0; i < R; i++){
		 out4<<w[i]<<"\t";
		}
		out4.close();

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

	//	int info = 0;
//		const char d = 'U';
//		const int N = 3;


//		lapack_int ch = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', 3, mat, 3);
	//	lapack_int inv = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', 3, mat, 3);
		 
		/*double amat[3][3] = {{0,1,3},{1,3,2},{3,2,-1}};
		double bmat[3][3] = {{5,-2,1},{-2,7,2},{1,2,6}};



		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				cout<<bmat[i][j]<<" ";
			}
		}
		
		

		lapack_int ch3 = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', 3, bmat[0], 3);
		
			for(int i = 0; i < 3; i++){
						for(int j = 0; j < 3; j++){
							cout<<bmat[i][j]<<" ";
						}
			}*/
		
		delete[] emat;
		delete[] gmat;
		delete[] w;
		delete[] bFunctions;

		}
	return 0;
}

int parity(int k, int l, int m, int i){
	double pvec[3] = {0,0,0};

	if(i == 0){
		pvec[0] =  pow(-1,k+1) ;
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


	if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == 1)){
	return 1;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
	return 2;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
	return 3;
	}
	else if( (pvec[0] == 1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
	return 4;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == 1)){
	return 5;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == 1)&&(pvec[2] == -1)){
	return 6;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == 1)){
	return 7;
	}
	else if( (pvec[0] == -1)&&(pvec[1] == -1)&&(pvec[2] == -1)){
	return 8;
	}
	else return -1;
}

int compPnumb(const void* b1, const void* b2){
	Basis::basisFunction * bas1 = (Basis::basisFunction*) b1;
	Basis::basisFunction * bas2 = (Basis::basisFunction*) b2;

	return 	(*bas1).pnumb - (*bas2).pnumb;

}

double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax){
	double intVal = 8;

	if( ( ( (*b1).xk + (*b2).xk)%2 == 1 )|| ( ( (*b1).yl + (*b2).yl)%2 == 1 ) || ( ( (*b1).zm + (*b2).zm)%2 == 1 ) ){
		intVal *= 0;
	}
	else{
		intVal *= ((1/( (double)(*b1).xk + (double)(*b2).xk + 1))*pow(xmax, (double)(*b1).xk + (double)(*b2).xk + 1))*((1/( (double)(*b1).yl + (double)(*b2).yl + 1))*pow(ymax, (double)(*b1).yl + (double)(*b2).yl + 1))*((1/( (double)(*b1).zm + (double)(*b2).zm + 1))*pow(zmax, (double)(*b1).zm + (double)(*b2).zm + 1));
	}
	return intVal;
}

double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax){
	double intVal = 8;

	double powers[3] = {(double)(*b1).xk + (double)(*b2).xk,   (double)(*b1).yl + (double)(*b2).yl,   (double)(*b1).zm + (double)(*b2).zm};
	
    int* ad1 = (int*) b1;
	int* ad2 = (int*) b2;

	ad1 += d1;
	ad2 += d2;

	intVal *= ((double) *ad1)*((double) *ad2);

	powers[d1]--;
	powers[d2]--;

	if((powers[0] <0 )||(powers[1] <0 )||(powers[2] <0 )){
		return 0;
	}

	if((( (int) powers[0]%2 == 1) ||( (int) powers[1]%2 == 1)||( (int) powers[2]%2 == 1)) ){
		return 0;
	}

	intVal *= (1/(powers[0]+1)*pow(xmax,powers[0]+1))*(1/(powers[1]+1)*pow(ymax,powers[1]+1))*(1/(powers[2]+1)*pow(zmax,powers[2]+1));
	
	return intVal;
}