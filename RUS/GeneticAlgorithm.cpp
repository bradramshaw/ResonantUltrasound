#include "stdafx.h"
#include "GeneticAlgorithm.h"

static const int nVars = 8; // number of variables, F, dF, etc...
static const int nParams = 9;  //total parameters including chiSq

static const int nThreads = 4;
double totalTime = 0;  // for timing things
LARGE_INTEGER freq;
_CrtMemState s1,s2,s3;

GeneticAlgorithm::GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability){	
	QueryPerformanceFrequency(&freq);
	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	initializeParameters(dataSet, dataSetLength, nPopulation, scaleFactor, crossingProbability);	
}

GeneticAlgorithm::~GeneticAlgorithm(void){
	delete [] _populationParametersNew;
	delete [] _populationParametersOld;
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
		_populationParametersOld[i].c11 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c22 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c33 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c44 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c55 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c66 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c12 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c13 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c23 = (randomDouble(0.0,1.0))*pow(10,9);
		
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
		_populationParametersOld[i].c22 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c33 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c44 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c55 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c66 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c12 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c13 = (randomDouble(0.0,1.0))*pow(10,9);
		_populationParametersOld[i].c23 = (randomDouble(0.0,1.0))*pow(10,9);
		
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
	std::cout<<"c11: "<<_minimumParameters.c11<<" "<<"c22: "<<_minimumParameters.c22<<" "<<"c33: "<<_minimumParameters.c33<<" "<<"c44: "<<_minimumParameters.c44<<" "<<"c55: "<<_minimumParameters.c55<<" "<<"c66: "<<_minimumParameters.c66<<" "<<std::endl<<"c12: "<<_minimumParameters.c12<<" "<<"c13: "<<_minimumParameters.c13<<" "<<"c23: "<<_minimumParameters.c23<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
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


double * GeneticAlgorithm::calculateFrequencies(double const * parameters){	
	double * frequencies;
	
	return frequencies;
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