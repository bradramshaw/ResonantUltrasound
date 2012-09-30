#pragma once
#include "Parameters.h"

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability);
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

	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double * _dataSet;
	
	double ** _residualArray;
	double ** _paramArray;

	Parameters::fitParameters * _populationParametersOld, *_populationParametersNew;
	Parameters::fitParameters _minimumParameters;

	void initializeRandomNumberGenerators();
	void initializeParameters(double* dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability);
	double randomDouble(double min, double max);
	double calculateResidual(Parameters::fitParameters * parameters,int threadID);
	double * calculateFrequencies(double const * parameters);
	
	void resetParameters(int nPopulation, double scaleFactor, double crossingProbability);
	static	UINT startResidualThread(LPVOID param);
	void residualCalculatingThread(Parameters::arrayBounds * arrayBounds);

	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm* pThis;
	};

};

