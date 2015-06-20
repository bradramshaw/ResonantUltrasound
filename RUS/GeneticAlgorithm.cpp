#include "stdafx.h"
#include "GeneticAlgorithm.h"




GeneticAlgorithm::GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int order, double xHL, double yHL, double zHL, double density, int nMissing){	
	
	totalTime = 0;  // for timing things
		
	QueryPerformanceFrequency(&freq);

	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	
	_nMissing = nMissing;
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
			for(int j = 0; j < _nMissing; j++){
				_minimumParameters.missFreq[j] = _populationParametersOld[i].missFreq[j];
			}
		}
	}
}



double GeneticAlgorithm::calculateResidual(Parameters::fitParameters * parameters, int threadID){
	
	double total = 0;	
	double * paramPointer = &(parameters->c11);	

	for(int i = 0; i < nParams; i++){
		_paramArray[threadID][i] = *paramPointer;
	
		paramPointer++;
	}

	
	double * frequencies = calculateFrequencies(_paramArray[threadID]);
		
	 

	if(_nMissing > 0){
	 std::vector<double> freqVect(0);


	 for(int i = 0; i < _dataSetLength; i ++){
		 freqVect.push_back(_dataSet[i]);
	 }
	 std::vector<double>::iterator iter;


	 double resTemp = INFINITY;
	 int totalLength = _nMissing + _dataSetLength;
	 long * missInd = new long[_nMissing];
	 bool * flag = new bool;
	 *flag = false;
	 long combs = nCombs(totalLength, _nMissing);

	 for(int i = 0; i < _nMissing; i++){
			missInd[i] = i;
		}
	
		for(int i = 0; i < combs; i ++){
			total = 0;
			std::vector<double> freqVectTemp(freqVect);
			iter = freqVectTemp.begin();
			
				for(int j = 0; j < _nMissing; j ++){
					iter += missInd[j];
					freqVectTemp.insert(iter, 0);
					iter = freqVectTemp.begin();
				}

				
				for(int m = 0; m < totalLength; m++){
					if(freqVectTemp[m] != 0){
			
					_residualArray[threadID][m] = (frequencies[m] - freqVectTemp[m])/(frequencies[m]);

						
					}
				
					else{
				
						_residualArray[threadID][m] = 0;
					}	
				
				}	
	
				for(int m = 0; m < totalLength; m++){
				
					total += _residualArray[threadID][m]*_residualArray[threadID][m];
				}
			
			

				if(total < resTemp){
					resTemp = total;
			
					for(int p =  0 ; p < _nMissing; p++){
						parameters->missFreq[p] = missInd[p];
					}
				}
			
				shiftInd(missInd, _nMissing -1, _nMissing, totalLength, flag);

			}
		
		total = resTemp;
	delete[] missInd;
	delete flag;
	}

	else{	
			for(int i = 0; i < _dataSetLength; i++){
				if(_dataSet[i]!= 0){
				_residualArray[threadID][i] = (frequencies[i] - _dataSet[i])/(frequencies[i]);
				}
				else{
					_residualArray[threadID][i] = 0;
				}
		
			}	
		
		for(int i = 0; i < _dataSetLength; i++){
			total += _residualArray[threadID][i]*_residualArray[threadID][i];
		}
	}


	delete[] frequencies;


	return total;
}

double * GeneticAlgorithm::calculateFrequencies(double * parameters){	

	
	
	double * frequencies = (double*) malloc(sizeof(double)*( _R - 6)) ;
	
	double **** ctens = initElasticConstants(parameters);

	

	double * emat = new double[_R*_R];
	memcpy(emat, _emat, sizeof(double)*_R*_R);

	
	double * gmat = calcGmat(_R, _basis, ctens, _gradientCalcs);

	double * temp = calcEigs(_R, emat, gmat);	
	
	
	for(int i  = 6; i < _R; i++){
		frequencies[i-6] = (sqrt(temp[i]))/(2*3.1415926535897);	
		
	}
	
	
	delete [] emat;
	delete [] gmat;
	delete [] temp;
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j++){
			for(int k = 0; k < 3; k++){
				delete[] ctens[i][j][k];
			}
			delete[] ctens[i][j];
		}
		delete[] ctens[i];
	}


	return frequencies;
}



void GeneticAlgorithm::calculateNewGenerations(int nGenerations){
	double averageTime = 0;
	

	for(int i = 0; i < nGenerations; i++){


	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints1, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints2, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints3, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, shuffleIndex, 0, _nPopulation );

	
	
	 //Parameters::fitParameters  temp;
	 //for(int r = 0; r < _nPopulation; r++){
		//temp = _populationParametersOld[r];
		//_populationParametersOld[r]  = _populationParametersOld[shuffleIndex[r]];
		//_populationParametersOld[shuffleIndex[r]] = temp;
	 //}
	
	 				
	for(int j = 0; j < _nPopulation; j++){

		
			double * pointerToOldVariable = &_populationParametersOld[j].c11;
			double * pointerToNewVariable = &_populationParametersNew[j].c11;
			double * pointerTog1Variable = &_populationParametersOld[ints1[j]].c11;
			double * pointerTog2Variable = &_populationParametersOld[ints2[j]].c11;
			double * pointerTog3Variable = &_populationParametersOld[ints3[j]].c11;	
				
			tetragonalParameters(pointerToOldVariable, pointerToNewVariable, pointerTog1Variable, pointerTog2Variable, pointerTog3Variable);

				//could be optimzed for vector arithmetic
			/*for(int k = 0; k < nVars; k++){
				double p = randomDouble(0,1);
				if(p > _crossingProbability){
					*pointerToNewVariable = *pointerToOldVariable;
				}
				else{
					*pointerToNewVariable = *pointerTog1Variable + _scaleFactor*(*pointerTog2Variable - *pointerTog3Variable);
				}
				pointerToNewVariable++;
				pointerToOldVariable++;
				pointerTog1Variable++;
				pointerTog2Variable++;
				pointerTog3Variable++;
			}	*/	
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
	
	//	std::cout<<s2.lTotalCount<<std::endl;

		WaitForMultipleObjects(nThreads,threadEvents,TRUE,INFINITE);	

		
		for(int timerIndex = 0; timerIndex < nThreads; timerIndex++){
		totalTime += threadContents[timerIndex].arrayBounds.time;
		}
		
		averageTime +=totalTime;
		
		mkl_free_buffers();
	
	/*	calculateMinimum();
		exportChiSq();*/
	}
	std::cout<<"Average time per generation for a thread: "<<averageTime/(nThreads*nGenerations)<<"ms"<<std::endl<<std::endl;
		
		
}

UINT GeneticAlgorithm::startResidualThread(LPVOID param){
	threadContents * contents = (threadContents*) param;
	contents->pThis->residualCalculatingThread(&(contents->arrayBounds));
	return 0;
}

void GeneticAlgorithm::residualCalculatingThread(Parameters::arrayBounds * arrayBounds){
	
		
	LARGE_INTEGER time1,time2;


	
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

//double ** GeneticAlgorithm::derivatives(){
//
//
//
//}

void GeneticAlgorithm::exportChiSq(){
    	std::ofstream out;
		out.open("chiSq3.dat",std::ios_base::app);
		out.precision(15);
		out<<_minimumParameters.chiSq<<"\t";
		out.close();
}

void GeneticAlgorithm::printMinimumParameters(){

		double rms = 100 * sqrt(_minimumParameters.chiSq/(_nMissing + _dataSetLength));

	std::cout<<std::left<<std::setfill(' ')<<std::setw(10)<<"c11"<<std::setw(10)<<"c22"<<std::setw(10)<<"c33"<<std::setw(10)<<"c13"<<std::setw(10)<<"c23"<<std::setw(10)<<"c12"<<std::setw(10)<<"c44"<<std::setw(10)<<"c55"<<std::setw(10)<<"c66"<<std::endl;
	std::cout<<std::left<<std::setfill(' ')<<std::setw(10)<<_minimumParameters.c11/(pow(10,9))<<std::setw(10)<<_minimumParameters.c22/(pow(10,9))<<std::setw(10)<<_minimumParameters.c33/(pow(10,9))<<std::setw(10)<<_minimumParameters.c13/(pow(10,9))<<std::setw(10)<<_minimumParameters.c23/(pow(10,9))<<std::setw(10)<<_minimumParameters.c12/(pow(10,9))<<std::setw(10)<<_minimumParameters.c44/(pow(10,9))<<std::setw(10)<<_minimumParameters.c55/(pow(10,9))<<std::setw(10)<<_minimumParameters.c66/(pow(10,9))<<std::endl<<std::endl<<"Residual: "<< rms <<" %"<<std::endl<<std::endl;

	//std::cout<<"c11: "<<_minimumParameters.c11<<" "<<"c22: "<<_minimumParameters.c22<<" "<<"c33: "<<_minimumParameters.c33<<std::endl<<"c44: "<<_minimumParameters.c44<<" "<<"c55: "<<_minimumParameters.c55<<" "<<"c66: "<<_minimumParameters.c66<<"c12: "<<_minimumParameters.c12<<" "<<"c13: "<<_minimumParameters.c13<<" "<<"c23: "<<_minimumParameters.c23<<std::endl<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
	
		std::ofstream out;
		out.open("output.dat");
		out.precision(15);
		out<<_minimumParameters.c11<<'\t'<<_minimumParameters.c22<<'\t'<<_minimumParameters.c33<<'\t'<<_minimumParameters.c44<<'\t'<<_minimumParameters.c55<<'\t'<<_minimumParameters.c66<<'\t'<<_minimumParameters.c12<<'\t'<<_minimumParameters.c13<<'\t'<<_minimumParameters.c23;
		out.close();


	double *frequencies = calculateFrequencies(&(_minimumParameters.c11));
		if(_nMissing > 0){
			std::vector<double> freqVect(0);
			 for(int i = 0; i < _dataSetLength; i ++){
			 freqVect.push_back(_dataSet[i]);
		 }
		 std::vector<double>::iterator iter = freqVect.begin();
		
			for(int j = 0; j < _nMissing; j ++){
				iter += _minimumParameters.missFreq[j];
				freqVect.insert(iter, 0);
				iter = freqVect.begin();
			}

			std::cout<<std::left<<std::setfill(' ')<<std::setw(20)<<"Measured"<<std::setw(20)<<"Calculated"<<std::setw(20)<<"%diff"<<std::endl;
	
			for(int i = 0; i < _dataSetLength+_nMissing; i ++){
				std::cout<<std::setw(20)<<freqVect[i]<<std::setw(20)<<frequencies[i]<<std::setw(20)<<(freqVect[i] - frequencies[i])/(frequencies[i]) * 100<<std::endl;
			}
			for(int i = _dataSetLength+_nMissing; i < _dataSetLength+_nMissing+ 10; i ++){
				std::cout<<std::setw(20)<<' '<<frequencies[i]<<std::endl;
		}
		}

		else{
		std::cout<<std::left<<std::setfill(' ')<<std::setw(20)<<"Measured"<<std::setw(20)<<"Calculated"<<std::setw(20)<<"%diff"<<std::endl;
		
		for(int i = 0; i < _dataSetLength; i ++){
			std::cout<<std::setw(20)<<_dataSet[i]<<std::setw(20)<<frequencies[i]<<std::setw(20)<<(_dataSet[i] - frequencies[i])/(frequencies[i]) * 100<<std::endl;
		}
		for(int i = _dataSetLength+_nMissing; i < _dataSetLength+_nMissing+ 10; i ++){
			std::cout<<std::setw(20)<<' '<<frequencies[i]<<std::endl;
		}

		}
	

		std::ofstream fout;
		fout.open("freqoutput.dat");
		fout.precision(15);
		for(int i = 0; i < _R; i ++){
			fout<<frequencies[i]<<std::endl;
		}
		fout.close();

		delete[] frequencies;

}