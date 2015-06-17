#pragma once
//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <sstream>
//#include <exception>

//herherhherherh

class DataExtractor
{

private:
	int dataReader(std::string name, std::vector<std::string>* buffer);	// Reads the data file into a vector in string format, returns number of lines
	bool dataConverter(std::vector<std::string>* stringDat, double* numDat);
	template<typename T> T* AllocateDynamicArray(int entries);
	template <typename T> void FreeDynamicArray(T* dArray);
	double* dataArray;
	int number_of_lines;

public:
	DataExtractor(std::string name);
	~DataExtractor();
	double* getDataArray();
	int getNumberOfLines();

};

