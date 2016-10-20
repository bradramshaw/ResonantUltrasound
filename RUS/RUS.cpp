// RUS.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
using namespace std;

// All of the non-zero elastic constants, in pascals. The notation is 1 xx, 2 yy, 3 zz, 4 yz, 5 xz, 6 xy. there are symmetries, and so for example c12 = c_xxyy = c_yyxx = c21

//
//static const double density = 7956; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.00205/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.00212/2; // HALF length in the y direction
//static const double zHL = 0.00025/2; // HALF length in the z direction. CeCoIn5****

//static const double density = 10079; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.00307/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.00278/2; // HALF length in the y direction
//static const double zHL = 0.00258/2; // HALF length in the z direction. URu2Si2****

//static const double density = 8174; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.00124/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.00212/2; // HALF length in the y direction
//static const double zHL = 0.00262/2; // HALF length in the z direction. Nb****

//static const double density = 6302; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.000856/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.001054/2; // HALF length in the y direction
//static const double zHL = 0.00018/2; // HALF length in the z direction. YBCO67****

//static const double density = 8422; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.0006384/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.0022416/2; // HALF length in the y direction
//static const double zHL = 0.0022154/2; // HALF length in the z direction. PuCoGa5***

//static const double density = 3501; // density of the material in grams/meter^3. All units are SI
//
//static const double xHL = 0.0015/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.0015/2; // HALF length in the y direction
//static const double zHL = 0.002/2; // HALF length in the z direction. Diamond****

static const double density = 6700; // density of the material in (Kg?) grams/meter^3. All units are SI

static const double xHL = 0.000823/2; // HALF length in the x direction, in meters.
static const double yHL = 0.000674/2; // HALF length in the y direction
static const double zHL = 0.000120/2; // HALF length in the z direction. Hg1201****

//static const double density = 6890; // density of the material in (Kg?) grams/meter^3. All units are SI
//
//static const double xHL = 0.001239/2; // HALF length in the x direction, in meters.
//static const double yHL = 0.00099/2; // HALF length in the y direction
//static const double zHL = 0.00033/2; // HALF length in the z direction. Hg1201****


	/*		QueryPerformanceCounter(&time1);*/
	
	//	QueryPerformanceCounter(&time2);
	//std::cout<<"Time to memcpy: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl;

int _tmain(int argc, _TCHAR* argv[]) //main function
{
	_CrtMemState s1,s2,s3;
	LARGE_INTEGER gtime1,gtime2,gfreq;  // stores times and CPU frequency for profiling
	QueryPerformanceFrequency(&gfreq);
	//double **** ctens = initElasticConstants(); // 4 dimensional elastic constant array. can probably be simpler (obviously)

	
		VSLStreamStatePtr stream;
		SYSTEMTIME t;
	    GetLocalTime(&t);
	    vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
	
		DataExtractor extractor("../RUS/HG1201_295K.dat");
		double * data = extractor.getDataArray();
		int nPoints = extractor.getNumberOfLines();
		
		int order, nMissing; // will store the max order of the polynomials to use
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

		cout << "Maximum number of missing peaks? ";
		cin >> nMissing;
		cout << endl;

		GeneticAlgorithm geneticAlgorithm(data, nPoints, 50, scale, cross, order, xHL, yHL, zHL, density, nMissing);
	
		geneticAlgorithm.calculateMinimum();
		geneticAlgorithm.printMinimumParameters();	

	while(true){ // bad programming...

		int nGens;
		cout<<"Number of generations: ";
		cin>>nGens;
		cout<<endl;

			//CURRENTLY ONE THREAD CHECK
		//_CrtMemCheckpoint( &s1 );

		QueryPerformanceCounter(&gtime1);
	
		geneticAlgorithm.calculateNewGenerations(nGens);	

		QueryPerformanceCounter(&gtime2);
		std::cout<<"Total time per generation: "<<(1000*(double)(gtime2.QuadPart-gtime1.QuadPart)/(gfreq.QuadPart))/nGens<<"ms"<<std::endl<<std::endl;

	/*	_CrtMemCheckpoint( &s2 );
		if ( _CrtMemDifference( &s3, &s1, &s2 ) )
      _CrtMemDumpStatistics( &s3 );

	cout<<s3.lTotalCount<<endl;*/


		geneticAlgorithm.calculateMinimum();	

		geneticAlgorithm.printMinimumParameters();	
		
		
	

		}
	return 0;
}
