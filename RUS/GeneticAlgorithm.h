#pragma once
#include "Parameters.h"
#include "Basis.h"

	
#define		timing1 QueryPerformanceCounter(&stime1)
#define		timing2 QueryPerformanceCounter(&stime2)
#define		timingOut	std::cout<<"Time is: "<<1000*(double)(stime2.QuadPart-stime1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int order, double xHL, double yHL, double zHL, double density, int nMissing);
	~GeneticAlgorithm(void);
	
	void calculateMinimum();
	void printMinimumParameters();
	void exportChiSq();
	void calculateNewGenerations(int nGenerations);
	
private:
	double totalTime;
	_CrtMemState s1,s2,s3;
	LARGE_INTEGER stime1,stime2,freq;  // stores times and CPU frequency for profiling

	static const int nVars = 9; // number of variables, F, dF, etc...
	static const int nParams = 10;  //total parameters including chiSq

	static const int nThreads = 8;

    VSLStreamStatePtr stream;
	int * ints1;
	int * ints2;
	int * ints3;
	int * shuffleIndex;

	int _order;
	int _R;
	Basis::basisFunction * _basis;

	double _xHL, _yHL, _zHL;
	double _density;
	int _nMissing;

	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double * _dataSet;

	double * _emat;

	int * _basisPop;
	
	double * _gradientCalcs;

	double ** _residualArray;
	double ** _paramArray;

	static int compPnumb(const void * b1, const void * b2); //comparitor for parity numbers. 
	Basis::basisFunction * createBasis(int order, int* basisPop);
	int parity(int k, int l, int m, int coord, int * basisPop); // parity function looks at the symmetry of a basis function, ie x^2 * y * z^3. more below in the full function definition.	
	double **** initElasticConstants(double* parameters); //initialize the full tensor. 

	void initialiseMatrices();
	double * calcEmat(int R,Basis::basisFunction * bFunctions);  // kinetic energy
	
	double * calcGmat(int R, Basis::basisFunction * bFunctions, double **** ctens, double * gradientCalcs); // elastic energy

	double * calcEigs(int R, double * emat, double * gmat);// eigenvalues (resonant frequencies squared, or maybe their inverse. anyway it's obvious)

	double * calcGradient(int R, Basis::basisFunction * bFunctions);
	double ** derivatives();

	void isotropicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);	
	void cubicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);
	void hexagonalParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);
	void tetragonalParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);
	void orthorhombicParameters(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);

	double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax); // This integrates  a pair of basis functions within the limits specified. Note that this assumes a parallelapiped, and takes "half" dimensions as inputs
	double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax); // integrates basis functions after differentiating with respect to one coordinate in each basis function of the pair. 

	static int dComp(const void *a, const void *b);


	Parameters::fitParameters * _populationParametersOld, *_populationParametersNew;
	Parameters::fitParameters _minimumParameters;

	void initializeRandomNumberGenerators();
	void initializeParameters(double* dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability);
	double randomDouble(double min, double max);
	double calculateResidual(Parameters::fitParameters * parameters,int threadID);
	double * calculateFrequencies(double * parameters);
	
	void resetParameters(int nPopulation, double scaleFactor, double crossingProbability);
	static	UINT startResidualThread(LPVOID param);
	void residualCalculatingThread(Parameters::arrayBounds * arrayBounds);
	

	long nCombs(int N, int k);
	void shiftInd(long * indV,long shiftIndex, long nMiss, long nFreq, bool * maxFlag);
	


	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm* pThis;
	};

};
