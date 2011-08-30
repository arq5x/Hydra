/********************************************************************************
Program: 		DnaLibrary.h
Description: 	Data structure for storing the details of a particular DNA library.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*********************************************************************************/
#ifndef DNALIBRARY_H
#define DNALIBRARY_H

using namespace std;

struct DNALIB {
    string _sample;
    string _mappingFile;
    float  _expSize;
    float  _variance;
    float  _numDevs;
    
    // constructor
	DNALIB (const string sample,
	        const string mappingFile,
	        const int size, 
	        const int variance, 
	        const int degreesOfVariance) :
		_sample(sample),
		_mappingFile(mappingFile),
		_expSize(size), 
		_variance(variance), 
		_numDevs(degreesOfVariance) 
	{}
    
    // get the length dev based on the description
    // of the library.
	float getMLD(void) const
	{ 
	    // should this be 2* below?
        return _numDevs * _variance;
	}

    // get the non-overlap based on the description
    // of the library.
	float getMNO(void) const
	{ 
        return 2 * (_expSize + (_numDevs * _variance));
	}	
};

#endif
