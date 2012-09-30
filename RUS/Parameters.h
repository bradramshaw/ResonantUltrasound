#pragma once
namespace Parameters
{
	__declspec( align( 128) ) struct fitParameters
	{
		double c11; //0
		double c22; //1
		double c33; //2
		double c44; //3
		double c55; //4
		double c66; //5
		double c12; //6
		double c13; //7
		double c23; //8
		double chiSq; //9
		
	};

    struct arrayBounds
	{
		double start;
		double end;
		HANDLE handle;
		double time;
		int threadID;
	};
}