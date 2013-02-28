#include <iostream>
#include <ctime>
//#include <omp.h>
#include "Hydra.h"
#include "Version.h"

using namespace std;

// define our program name
const string PROGRAM_NAME  = "hydra";

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
bool IsNumber(const std::string& s);

int main(int argc, char* argv[]) {

    int memory = (int) ((8) * pow(2,30));  // 1Gb = 2^30, use 8Gb by default.
    int minSupport        = 2;
    int maxLinkedDistance = 1000000; // 1Mb by default
    int numThreads = 1;
    int maxMappings = 100000;
    
    vector<DNALIB> sampleLibs;

	// sample ids
    vector<string> samples;    
    // input files
    vector<string> discordantFiles;
    
    // library parameters for each file.
    vector<int> expectedSizes;
    vector<int> variances;
    vector<int> degreesOfVariance;

    bool haveLD                 = false;  // deprecated 
    bool haveSD                 = false;  // deprecated
    bool haveSizes              = false;
    bool haveVariances          = false;
    bool haveDegressOfVariances = false;
    bool haveMinSupport         = false;
    bool haveMaxLinkedDistance  = false;
    bool lumpInversions         = false;
    bool ignoreSize             = false;

    string mappingUsage         = "best";
    unsigned int editBeyondBest = 0;


    // initialization
    ostringstream sb;
    char* end_ptr = NULL;

    // our configuration variables
    bool showHelp          = false;

    // output files
    string configFile;
    string routedFile;
    
    // checks for existence of parameters
    bool haveConfigFile         = false;
    bool haveRoutedFile         = false;
    

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-h", 2, parameterLength) || 
        PARAMETER_CHECK("--help", 5, parameterLength)) {
            showHelp = true;
        }
    }

    if(showHelp) ShowHelp();

    // do some parsing (all of these parameters require 2 strings)
    cerr << endl << "Parameters:" << endl;
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-config", 7, parameterLength)) {
            if ((i+1) < argc) {
                haveConfigFile = true;
                configFile = argv[i + 1];
                cerr << "  Configuration file (-config): " << configFile << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-routed", 7, parameterLength)) {
            if ((i+1) < argc) {
                haveRoutedFile = true;
                routedFile     = argv[i + 1];
                cerr << "  Using routed file as input: " << routedFile << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-maxMappings", 12, parameterLength)) {
            if ((i+1) < argc) {
                maxMappings     = atoi(argv[i + 1]);
                cerr << "  Maximum mappings allowed before \"punting\": " << maxMappings << endl;
                i++;
            }
        }
        else {
            cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }       
    }

    if (!haveConfigFile) {
        cerr << "*****ERROR: You must specify an input configuration file.*****" << endl << endl;
        showHelp = true;
    }
    if (!haveRoutedFile) {
        cerr << "*****ERROR: You must specify an input routed file.*****" << endl << endl;
        showHelp = true;
    }
    
    if (!showHelp) {

        cerr << endl << "Processing: " << endl;
        
        if (haveConfigFile == true) {
            // open the mapping files for reading
            ifstream config(configFile.c_str(), ios::in);

            if ( !config ) {
                cerr << "Error: The configuration file (" << configFile << ") could not be opened. Exiting!" << endl;
                exit (1);
            }
            else {
                string sample, file;
                float  mean, std, numstd;
                while (config >> sample >> file >> mean >> std >> numstd) {
                    DNALIB lib(sample, file, mean, std, numstd);
                    sampleLibs.push_back(lib);
                }
            }
        }
        
        HydraPE *events = new HydraPE(sampleLibs, "", minSupport,
                                                maxLinkedDistance, ignoreSize, lumpInversions,
                                                mappingUsage, editBeyondBest, memory, true);
        
        events->LoadRoutedFile(routedFile);
        events->SortAllMasterFilesByPosition_New();
        events->FindPositionClusters_New();
        events->AssembleClusters(maxMappings);

        return 0;
    }
    else {
        ShowHelp();
    }
}


void ShowHelp(void) {
    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
    
    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Calls SV breakpoints from discordant paired-end mappings." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -config FILE -routed FILE -out STRING" << endl << endl;

    cerr << "Options:" << endl;
    cerr << "  -config\tConfiguration file." << endl;
    cerr << "       \t\tCol 1. Sample Id (string)" << endl;
    cerr << "       \t\tCol 2. Mapping file (path/file)" << endl;
    cerr << "       \t\tCol 3. Expected insert size (integer)" << endl;
    cerr << "       \t\tCol 4. Variance (integer)" << endl;
    cerr << "       \t\tCol 5. Num. variances (integer)" << endl << endl;
    
    cerr << "  -routed\tA single routed chr/chr/strand/strand file from HydraRouter." << endl << endl;
        
    cerr << "  -maxMappings\tMaximum number of mappings in a cluster before Hydra will \"punt\".." << endl << endl;
        
    // end the program here
    exit(1);
}


bool IsNumber(const std::string& s) {
   for (int i = 0; i < s.length(); i++) {
       if (!std::isdigit(s[i]))
           return false;
   }
   return true;
}

