#pragma once
#include "Parameters.h"
#include "Basis.h"

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int order, double xHL, double yHL, double zHL, double density);
	~GeneticAlgorithm(void);
	
	void calculateMinimum();
	void printMinimumParameters();
	void exportChiSq();
	void calculateNewGenerations(int nGenerations);
	
private:

    VSLStreamStatePtr stream;
	int * ints1;
	int * ints2;
	int * ints3;

	int _order;
	int _R;
	Basis::basisFunction * _basis;

	double _xHL, _yHL, _zHL;
	double _density;

	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double * _dataSet;
	
	double ** _residualArray;
	double ** _paramArray;

	static int compPnumb(const void * b1, const void * b2); //comparitor for parity numbers. 
	Basis::basisFunction * createBasis(int order);
	int parity(int k, int l, int m, int coord); // parity function looks at the symmetry of a basis function, ie x^2 * y * z^3. more below in the full function definition.	
	double **** initElasticConstants(double* parameters); //initialize the full tensor. 

	double * calcEmat(int R,Basis::basisFunction * bFunctions);  // kinetic energy
	double * calcGmat(int R, Basis::basisFunction * bFunctions, double **** ctens); // elastic energy
	double * calcEigs(int R, double * emat, double * gmat);// eigenvalues (resonant frequencies squared, or maybe their inverse. anyway it's obvious)

	double integrateBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, double xmax, double ymax, double zmax); // This integrates  a pair of basis functions within the limits specified. Note that this assumes a parallelapiped, and takes "half" dimensions as inputs
	double integrateGradBasis(Basis::basisFunction * b1, Basis::basisFunction * b2, int d1, int d2, double xmax, double ymax, double zmax); // integrates basis functions after differentiating with respect to one coordinate in each basis function of the pair. 



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

	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm* pThis;
	};

};

