#pragma once


class LoadData
{

private:
	int dataReader(std::string name, std::vector<std::string>* buffer);	// Reads the data file into a vector in string format, returns number of lines
	bool dataConverter(std::vector<std::string>* stringDat, double* numDat);
	template<typename T> T* AllocateDynamicVector(int entries);
	template<typename T> void FreeDynamicVector(T* dArray);
	double* dataVector;
	int number_of_lines;

public:
	LoadData(std::string name);
	~LoadData();
	double* getDataVector();
	int getNumberOfLines();
};
