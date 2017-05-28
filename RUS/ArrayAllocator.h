#pragma once

namespace ArrayAllocator
{
	template<typename T>  T*  AllocateDynamicArray(int entries) //allocates an array
	{
		 T *dynamicArray;

		  dynamicArray = new T[entries];
		
		  return dynamicArray;
	} //allocates an array
	
	 template <typename T>  void  FreeDynamicArray(T* dArray) // deallocates the same
	{
		
		delete [] dArray;
	}
};
