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

    // initialization
    ostringstream sb;
    char* end_ptr = NULL;

    // our configuration variables
    bool showHelp          = false;

    int memory = (int) ((8) * pow(2,30));  // 1Gb = 2^30, use 8Gb by default.
    int minSupport        = 2;
    int maxLinkedDistance = 1000000; // 1Mb by default
    int numThreads = 1;
    
    // TAB-separated configuration file listing the samples, files and their statistics.
    // 1. sample ID (string)
    // 2. path/filename (string)
    // 3. median/mean (float)
    // 4. mad/std (float)
    // 5. num mads/stds. (float)
    string configFile;
    string routedFiles;
    string posSortedFiles;
    vector<DNALIB> sampleLibs;

    // sample ids
    vector<string> samples;    
    // input files
    vector<string> discordantFiles;
    
    // library parameters for each file.
    vector<int> expectedSizes;
    vector<int> variances;
    vector<int> degreesOfVariance;

    // output files
    string outFile;

    // checks for existence of parameters
    bool haveConfigFile         = false;    
    bool haveDiscordants        = false;
    bool haveOutFile            = false;
    bool haveRoutedFiles        = false;
    bool havePosSortedFiles    = false;
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
        else if(PARAMETER_CHECK("-mem", 4, parameterLength)) {
            if ((i+1) < argc) {
                cerr << "  Memory (-mem): " << atoi(argv[i + 1]) << " Gb." << endl;
                // convert to bytes.
                memory = (int) (atoi(argv[i + 1]) * pow(2,30));
                i++;
            }
        }
        else if(PARAMETER_CHECK("-t", 2, parameterLength)) {
            if ((i+1) < argc) {
                cerr << "  Num. threads (-t): " << atoi(argv[i + 1]) << endl;
                // convert to bytes.
                numThreads = atoi(argv[i + 1]);
                //omp_set_num_threads(numThreads);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-out", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveOutFile = true;
                outFile = argv[i + 1];
                cerr << "  Summary breakpoint output file (-out): " << outFile << endl;
                cerr << "  Detailed breakpoint output file: " << outFile << ".detail" << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-routedFiles", 12, parameterLength)) {
            if ((i+1) < argc) {
                haveRoutedFiles = true;
                routedFiles = argv[i + 1];
                cerr << "  Using already-routed files list in: " << routedFiles << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-posSortedFiles", 15, parameterLength)) {
            if ((i+1) < argc) {
                havePosSortedFiles = true;
                posSortedFiles = argv[i + 1];
                cerr << "  Using position clusters files list in: " << posSortedFiles << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ms", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveMinSupport = true;
                minSupport = atoi(argv[i + 1]);
                cerr << "  Minimum number of supporting pairs (-ms): " << minSupport << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-lnk", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveMaxLinkedDistance = true;
                maxLinkedDistance = atoi(argv[i + 1]);
                cerr << "  Maximum intrachromosomal distance before unlinked (-maxIntra): " << maxLinkedDistance << endl;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-is", 3, parameterLength)) {
            ignoreSize = true;
            cerr << "  Break cluster ties based on edit distance instead of size. " << endl;
        }
        else if(PARAMETER_CHECK("-li", 3, parameterLength)) {
            lumpInversions = true;
            cerr << "  Combine +/+ and -/- mappings when finding inversions. " << endl;
        }
        else if(PARAMETER_CHECK("-use", 4, parameterLength)) {
            
            if ((i+1) < argc) {
                string useValue = argv[i + 1];
                if (useValue == "best") {
                    mappingUsage = "best";
                    editBeyondBest = 0;
                    cerr << "  Using *best* mappings." << endl;
                }
                else if (useValue == "all") {
                    mappingUsage = "all";
                    editBeyondBest = INT_MAX;
                    cerr << "  Using *all* mappings." << endl;
                }
                else if ( IsNumber(useValue) ) {
                    mappingUsage = "withinBest";
                    editBeyondBest = atoi(useValue.c_str());
                    cerr << "  Using all mappings within " << editBeyondBest << "edits of best." << endl;                   
                }
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
    if (!haveOutFile) {
        cerr << "*****ERROR: You must specify an output file.*****" << endl << endl;
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
        
        HydraPE *events = new HydraPE(sampleLibs, minSupport,
                                                maxLinkedDistance, ignoreSize, lumpInversions,
                                                mappingUsage, editBeyondBest, memory);
        
        cerr << "  Routing discordant mappings to master chrom/chrom/strand/strand files." << endl;  
        events->RouteDiscordantMappings();
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

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -in FILE1 FILE2 ... FILEn -out <breakpoints> -mld <bp> -mno <bp>" << endl << endl;

    cerr << "Options:" << endl;
    cerr << "  -config\tConfiguration file." << endl;
    cerr << "       \t\tCol 1. Sample Id (string)" << endl;
    cerr << "       \t\tCol 2. Mapping file (path/file)" << endl;
    cerr << "       \t\tCol 3. Expected insert size (integer)" << endl;
    cerr << "       \t\tCol 4. Variance (integer)" << endl;
    cerr << "       \t\tCol 5. Num. variances (integer)" << endl << endl;

    cerr << "  -out\tStub for output files." << endl << endl;

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

