/*
  ***************************************************************************
   TestMerge.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.
 ***************************************************************************
*/

#define PROGRAM_NAME "mergeSort"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

// local includes
#include "Hydra.h"
#include "Sort.h"
#include "Ancillary.h"
#include "ExternalMergeSort.h"


// function declarations
void ShowHelp(void);
	

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp      = false;
	bool haveInFile    = false;
    bool haveOutFile   = false;
    bool haveChunkSize = false;
    bool haveTempPath  = false;
    bool compress      = false;	
	
    
	// input and output files
	string inFile, outFile;
    string tempPath = "";
    float chunkSize = 1;  // default to 1Gb of RAM

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();


	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

 		if (PARAMETER_CHECK("-i", 2, parameterLength)) {
			haveInFile = true;
			inFile = argv[i + 1];
			i++;
		}
		else if (PARAMETER_CHECK("-o", 2, parameterLength)) {
			haveOutFile = true;
			outFile = argv[i + 1];
			i++;
		}
		else if (PARAMETER_CHECK("-zip", 4, parameterLength)) {
			compress = true;
		}
		else if (PARAMETER_CHECK("-mem", 4, parameterLength)) {
			haveChunkSize = true;
			chunkSize = atof(argv[i+1]);
			i++;
        }
        else if (PARAMETER_CHECK("-tmp", 4, parameterLength)) {
            haveTempPath = true;
            tempPath = argv[i+1];
            i++;
        }
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}
	if (!haveInFile) {
		cerr << "*****ERROR: You must specify an input file (-i).*****" << endl << endl;
		showHelp = true;
	}
	
	if (!showHelp) {
        int chunkSizeInBytes = (int) ((chunkSize) * pow(2,30));  // 1Gb = 2^30
	    ExtMergeSort<PAIR> *merge = new ExtMergeSort<PAIR> (inFile, byStart1Desc, chunkSizeInBytes, compress, tempPath);
        merge->DivideAndSort();
        if (haveOutFile) {
            if (compress == false) {
                ofstream output(outFile.c_str(), ios::out);
                merge->Merge(&output);
            }
            else {
                ogzstream output(outFile.c_str(), ios::out);
                merge->Merge(&output);
            }
        }
        else {
            merge->Merge(&cout);
        }            
		delete merge;
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {
	
	cerr << endl << "Program: " << PROGRAM_NAME << " (v1.0)" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;
	
	cerr << "Summary: k-way, memory-assisted merge sort." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i FILE" << endl << endl;
	
	cerr << "Options: " << endl;
	cerr << "\t-o FILE\t\t"            	<< "Name of an output file.  Default is to write to stdout." << endl << endl;

	cerr << "\t-zip BOOL\t"            	<< "GZIP compress the temp and output files." << endl << endl;

	cerr << "\t-mem FLOAT\t"            << "Number of gigabytes to use for sorting. Default = 1Gb" << endl << endl;
	
	exit(1);	
}