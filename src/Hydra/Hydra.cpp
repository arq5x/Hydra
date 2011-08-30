/******************************************************************************
Program:        hydra pem
Description:    A comprehensize algorithm for detecting structural 
                variation in all genomic regions using paired-end mapping.
                    
Author:         Aaron Quinlan, Ph.D
                University of Virginia
                aaronquinlan@gmail.com
*******************************************************************************/
#include "Hydra.h"
#include "Sort.h"
#include <ctime>
#include "Ancillary.h"
#include "ExternalMergeSort.h"

// make
HydraPE::HydraPE(vector<DNALIB> libraries, 
                 int minSupport, int maxLinkedDistance, bool ignoreSize, 
                 bool lumpInversions, string mappingUsage, int editBeyondBest, int memory) 
: _libraries(libraries)
, minSupport(minSupport)
, maxLinkedDistance(maxLinkedDistance)
, ignoreSize(ignoreSize)
, lumpInversions(lumpInversions)
, mappingUsage(mappingUsage)
, editBeyondBest(editBeyondBest)
, _memory(memory)
{}


// and break
HydraPE::~HydraPE(void) {}


int HydraPE::GetLengthDeviation(const PAIR &pair) {
    return  2 * (_libraries[pair.fileNum]._numDevs 
                * 
                 _libraries[pair.fileNum]._variance);
}


int HydraPE::GetSpanDeviation(const PAIR &pair) {
    return 2 * (_libraries[pair.fileNum]._expSize 
               + 
               (_libraries[pair.fileNum]._numDevs * _libraries[pair.fileNum]._variance));
}


void HydraPE::WriteFragClusterToFile(ofstream *out, vector<PAIR> &cluster, int clusterId) {

    for (size_t i = 0; i < cluster.size(); ++i) {
        cluster[i].clusterId = clusterId;
        *out << cluster[i];
    }

}


void HydraPE::WritePosClusterToFile(ofstream *out, vector<PAIR> &cluster, int clusterId) {

    for (size_t i = 0; i < cluster.size(); ++i) {
        cluster[i].clusterId = clusterId;
        *out << cluster[i];
    }

}


/******************************************************************************
* RouteDiscordantMappings
*
* SUMMARY:
******************************************************************************/

void HydraPE::RouteDiscordantMappings() {
    
    // first check if we can open all of the requested files
    for (int fileNum = 0; fileNum < _libraries.size(); ++fileNum) {
        // open the mapping files for reading
        ifstream mappings(_libraries[fileNum]._mappingFile.c_str(), ios::in);        
        if ( !mappings ) {
            cerr << "Error: The requested mappings file (" << _libraries[fileNum]._mappingFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            cerr << "Found " << _libraries[fileNum]._mappingFile << endl;
            mappings.close();
        }
    }

    int numFiles = _libraries.size();
    string file; 
    for (int fileNum = 0; fileNum < numFiles; ++fileNum) {      
        file = _libraries[fileNum]._mappingFile;
        RouteFile(file, fileNum);
    }  // end loop through discordant mapping files.
}


void HydraPE::RouteFile(const string &file, int fileNum) {
    
    string currRead = "";
    string prevRead = "";
    CORE_PAIR line;
    PAIR      fullPair;

    // vector of mappings for the current read.
    pairVector currReadMappings;
    vector< pairVector> buffer;
    // reserve some memory so we don't have to reallocate (as much)
    time_t begin,end;
    // open the mapping files for reading
    ifstream *mappings;
    currReadMappings.reserve(1000);
    buffer.reserve(BUFFER_SIZE);

    time(&begin);       
    cerr << "Routing mappings from: " << file << "...";

    // open the mapping files for reading
    mappings = new ifstream(file.c_str(), ios::in);
     
    // loop through the current mapping file and store each line in a PAIR struct.
    while (*mappings >> line) {
        // get the id for the current mapping                           
        currRead = line.readId; 
        // set the ends of the pair
        // ****NOTE****: This assumes that whichMateIsBlock1 is in the Hydra input.
        if (line.whichMateIsBlock1 == 1) {
            line.mate1 = 1; line.mate2 = 2;
        }
        else {
            line.mate1 = 2; line.mate2 = 1;             
        }    
        
        fullPair.core = line;
        // initialize the other fields for the PAIR structure
        fullPair.support              = 0;
        fullPair.clusterId            = -1;   // initially, no mapping is part of a cluster.
        fullPair.used                 = false;
        fullPair.include              = true;
        fullPair.fileNum              = fileNum;

        // ensure that the mappings for each end of the pair are in the correct order.
        CorrectMateOrder(fullPair.core);

        // compute the diff b/w the mapping distance and the exp. size for this library.
        fullPair.diffFromExpectedSize = abs(getFragSize(fullPair.core) - (int) _libraries[fileNum]._expSize);

        // Add this mapping to the list of mappings for this READ.
        // We will use all of the mappings for a read to reduce the number of 
        // mappings to those that have the least number of mismatches.    
        if ( (currRead != prevRead) && (prevRead != "") ) {
        
            // reduce the mappings for this pair to those with the least
            // edit distance and add the remaning mappings to
            CullMappingsByMisMatches(currReadMappings);
            buffer.push_back(currReadMappings);
            
            if (buffer.size() == BUFFER_SIZE) {
                AddMappingsToMasterFile(buffer);
                buffer.clear();
            }

            // reset 
            currReadMappings.clear();
            currReadMappings.push_back(fullPair);
        }
        else {
            currReadMappings.push_back(fullPair);
        }
        prevRead = currRead;
    }
    // add the best mappings for the last pair in the file.
    CullMappingsByMisMatches(currReadMappings);
    // add the mappings to the appropriate chrom/chrom/strand/strand file
    buffer.push_back(currReadMappings);
    AddMappingsToMasterFile(buffer);
    buffer.clear();
    delete mappings;

    time(&end);
    double dif = difftime (end,begin);
    
    cout << "Time elapsed: " << dif << " sec"<< endl;   
}


void HydraPE::LoadRoutedFile(const string &routedFile) {
    // NEED LOGIC TO DETERMINE IF FILE IS INTER OR INTRA BASED ON THE FILENAME
    bool intra = false;
    int end_of_chroms = routedFile.size() - 2;

    size_t first_chrom_idx, second_chrom_idx;
    string first_chrom, second_chrom;
    first_chrom_idx  = routedFile.find("chr");
    second_chrom_idx = routedFile.find("chr", first_chrom_idx+1);
    first_chrom  = routedFile.substr(3, second_chrom_idx-3);
    second_chrom = routedFile.substr(second_chrom_idx+3, end_of_chroms-(second_chrom_idx+3));
    
    if (first_chrom == second_chrom) {
        intra = true;
    }
    _masterChromStrandFiles[routedFile]           = 0;
    _masterChromStrandFileTypes[routedFile]       = intra; 
}

void HydraPE::LoadRoutedFileList(const string &routedFiles) {
    ifstream routed(routedFiles.c_str(), ios::in);
    string file;
    bool intra;
    while (routed >> file >> intra ) 
    {
        _masterChromStrandFiles[file]     = 0;
        _masterChromStrandFileTypes[file] = intra; 
    }
}

void HydraPE::LoadPosSortedFileList(const string &sortedFiles) {
    ifstream sorted(sortedFiles.c_str(), ios::in);
    string file;
    while (sorted >> file)
        _posSortFiles.push_back(file);
}

void HydraPE::LoadPosClusterFileList(const string &posClusterFiles) {
    ifstream clusters(posClusterFiles.c_str(), ios::in);
    string file;
    while (clusters >> file)
        _posClusterFiles.push_back(file);
}


/******************************************************************************
* CorrectMateOrder
*
* Swaps the ends of a pair so that:
*    1. the leftmost end is first when an intra-chromosomal.
*    2. the lexigraphically lower chrom is first when an inter-chromosomal.
******************************************************************************/
void HydraPE::CorrectMateOrder(CORE_PAIR &pair) {
    if (pair.chrom1 == pair.chrom2) {
        if (pair.start2 < pair.start1) {
            SwapEnds(pair);
        }
    }
    else {
        if (pair.chrom2 < pair.chrom1) {
            SwapEnds(pair);
        }   
    }
}


/******************************************************************************
* SwapEnds
*
* Swaps the ends of a pair.
******************************************************************************/
void HydraPE::SwapEnds(CORE_PAIR &pair) {

    // hold end1 in temp
    string tmpChrom        = pair.chrom1;
    short tmpMate          = pair.mate1;
    int tmpStart           = pair.start1;
    int tmpEnd             = pair.end1;
    string tmpStrand       = pair.strand1;
    unsigned short tmpEdit = pair.edit1;
    
    // set end1 to end2
    pair.chrom1  = pair.chrom2;
    pair.mate1   = pair.mate2;
    pair.start1  = pair.start2;
    pair.end1    = pair.end2;
    pair.strand1 = pair.strand2;
    pair.edit1   = pair.edit2;
    
    // set end2 to end1
    pair.chrom2  = tmpChrom;
    pair.mate2   = tmpMate;
    pair.start2  = tmpStart;
    pair.end2    = tmpEnd;
    pair.strand2 = tmpStrand;
    pair.edit2   = tmpEdit;
}


/******************************************************************************
* AddMappingsToMasterMap
*
* Creates a map where the key is chrom1,chrom2,strand1,strand2 and the value
* is a vector of all mappings with the same chroms and strands.  In effect, 
* constructing this map does much of the sorting work for us ahead of time.
******************************************************************************/
void HydraPE::AddMappingsToMasterFile(const vector<pairVector> &readMappings) {

    vector<pairVector>::const_iterator readIter = readMappings.begin();
    vector<pairVector>::const_iterator readEnd  = readMappings.end();
    for (; readIter != readEnd; ++readIter) {
        
        pairVector::const_iterator mapIter = (*readIter).begin();
        pairVector::const_iterator mapEnd  = (*readIter).end();
        // write the contents of the current buffer to the master file
        for (; mapIter != mapEnd; ++mapIter) {
            // determine the name of the master file to write to.
            stringstream masterFileNameSS;
            masterFileNameSS << mapIter->core.chrom1 << mapIter->core.chrom2;
            if (this->lumpInversions == true && IsInversionMapping(mapIter->core) == true)
                masterFileNameSS << "+,+";
            else
                masterFileNameSS << mapIter->core.strand1 << mapIter->core.strand2;
            string masterFileName = masterFileNameSS.str();
        

            // setup the appropriate file handle to write the data to it
            bool isIntrachromosomal = true;
            if (mapIter->core.chrom1 != mapIter->core.chrom2)
                isIntrachromosomal = false;
                
            ostream *output = GetMasterFileStream(masterFileName, isIntrachromosomal);
            
            // require that the mapping has an inner span.
            // this prevents, for example, 150bp fragments with 100bp reads from being considered.
            // mol. .......................................
            // end1 11111111111111111111
            // end2                 22222222222222222222222
            // 
            // there is no span here, so w/o split-read mapping, we can't call a breakpoint.
            if (mapIter->core.start2 > mapIter->core.end1)
                *output << *mapIter;
        }
    }
}


ostream* HydraPE::GetMasterFileStream(const string &masterFileName, bool isIntra) {
    ostream *output;
    // we've already opened this master file.  just get it's file pointer.
    if (_masterChromStrandFiles.find(masterFileName) != _masterChromStrandFiles.end()) {
        output = _masterChromStrandFiles[masterFileName];
    }
    // we've not yet opened a file for this chr/chr/str/str combo.  do so and store the pointer.
    else {
        output = new ofstream(masterFileName.c_str(), ios::out);
        _masterChromStrandFiles[masterFileName]     = output;
        _masterChromStrandFileTypes[masterFileName] = isIntra;
    }
    return output;
}


bool HydraPE::IsInversionMapping(const CORE_PAIR &pair) {
    if (
        (pair.strand1 == "+" && pair.strand2 == "+") 
        || 
        (pair.strand1 == "-" && pair.strand2 == "-")
       ) 
    {
           return true;
    }
    return false;
}


void HydraPE::RemoveMasterChromStrandFiles() {
    cerr << "  Cleaning up old files." << endl;           
    
    map<string, ostream*>::const_iterator masterFileItr = _masterChromStrandFiles.begin();
    map<string, ostream*>::const_iterator masterFileEnd = _masterChromStrandFiles.end();
    for (; masterFileItr != masterFileEnd; ++masterFileItr) {
        string     file = masterFileItr->first;
        ostream *stream = masterFileItr->second;
        delete stream;
        remove(file.c_str());  // remove = UNIX "rm"
    }
}


void HydraPE::RemovePositionSortedFiles() {
    cerr << "  Cleaning up old files." << endl;           
    
    vector<string>::const_iterator fragFileItr = _posSortFiles.begin();
    vector<string>::const_iterator fragFileEnd = _posSortFiles.end();
    for (; fragFileItr != fragFileEnd; ++fragFileItr) {
        remove(fragFileItr->c_str());  // remove = UNIX "rm"
    }
}


void HydraPE::RemovePositionClusterFiles() {
    cerr << "  Cleaning up old files." << endl;           
    
    vector<string>::const_iterator posFileItr = _posClusterFiles.begin();
    vector<string>::const_iterator posFileEnd = _posClusterFiles.end();
    for (; posFileItr != posFileEnd; ++posFileItr) {
        remove(posFileItr->c_str());  // remove = UNIX "rm"
    }
}



void HydraPE::SortAllMasterFilesByPosition_New() {

    cerr << "Sorting groups by position." << endl;           
    
    // list of master files to sort
    vector<string> masterFiles;

    // extract the file names into our list.
    // NOTE: we are loading into a new vector so that
    // OpenMP can work on a for loop with an index, 
    // rather than an STL iterator (the latter requires >= GCC 4.4)
    map<string, ostream*>::const_iterator mItr = _masterChromStrandFiles.begin();
    map<string, ostream*>::const_iterator mEnd = _masterChromStrandFiles.end();
    for( ; mItr != mEnd; ++mItr ) 
    {
        masterFiles.push_back( mItr->first );
        bool is_intra = _masterChromStrandFileTypes[mItr->first];
        _posSortFileTypes.push_back(is_intra);
    }

    
    // loop through each master chrom/chrom/str/str file and sort it
    // by the differences from the expected frag. size.  name, compress 
    // and write the sorted file.
    //#pragma omp parallel for
    for (int i = 0; i < masterFiles.size(); ++i) 
    {
        string file = masterFiles[i];
     
        // log what we're doing
        time_t begin,end;
        time(&begin);
        cerr << "\tSorting " << file << " by position...";
        // create a sorted output file for the master file.
        string outSorted = file;
        outSorted.append(".posSorted");
        ofstream output(outSorted.c_str(), ios::out);

        // create new merge sorting object
        ExtMergeSort<PAIR> *merge = new ExtMergeSort<PAIR> (file, &output, byStart1Asc, _memory);

        // sort the master input file
        merge->DivideAndSort();

        // store the name of the sorted file
        _posSortFiles.push_back(outSorted);

        // write the sorted input to the output file.
        merge->Merge();
        output.close();
        delete merge;  
        
        time(&end);
        double dif = difftime (end,begin);
        
        cout << "Time elapsed: " << dif << " sec"<< endl;
    }
}


void HydraPE::FindPositionClusters_New() {

    cerr << "Finding possible breakpoint clusters by position." << endl;
    
    // loop through each chrom/chrom/str/str file that has been
    // sorted by fragments size and search for clusters of 
    // discordant mappings that have the same size within the
    // tolerances of their respective libaries.

    //#pragma omp parallel for
    for (int i = 0; i < _posSortFiles.size(); ++i) {

        // open the mapping file for reading
        string sortedFile  = _posSortFiles[i];
        ifstream *sortedMappings = new ifstream(sortedFile.c_str(), ios::in);

        // log what we're doing
        cerr << "\tFinding potential clusters in  " << sortedFile << "..." << endl;
        time_t begin,end;
        time(&begin);
         
        // create a file for storing potential clusters.  store it's name for later.
        string clusterFile = _posSortFiles[i];
        clusterFile.append(".posClusters");
        _posClusterFiles.push_back(clusterFile);
        ofstream *clusters = new ofstream(clusterFile.c_str(), ios::out);

        // look for clusters in the current input file.  
        int clusterId = 0;
        if (_posSortFileTypes[i] == true) // true means intra-chromosomal
        {
            FindIntraClusters(sortedMappings, clusters, clusterId);
        }
        else 
        {
            FindInterClusters(sortedMappings, clusters, clusterId);
        }
        
        // close the inout and output files    
        clusters->close();
        sortedMappings->close();
        delete clusters;
        delete sortedMappings;
        
        time(&end);
        double dif = difftime (end,begin);
        cout << "Time elapsed: " << dif << " sec"<< endl;
    }
}


void HydraPE::FindIntraClusters(ifstream *sortedMappings, ofstream *clusters, int &clusterId) 
{
    PAIR prevLine;
    PAIR currLine;
    vector<PAIR> currCluster;
    
    *sortedMappings >> prevLine;
    currCluster.push_back(prevLine);
    while (*sortedMappings >> currLine) 
    {
        /*
             p []---------[]
             c          []---------[]
        */
    
        // NEED separate logic for INTER and INTRA
        // while ((currLine.core.end1 < prevLine.core.start2) && 
        //        ((currLine.core.start1 - prevLine.core.end1) <= 2000)
        //       ) 
        while ((currLine.core.end1 < prevLine.core.start2) &&
               (sortedMappings->good())
              )             
        {
            currCluster.push_back(currLine);
            prevLine = currLine;
            *sortedMappings >> currLine;
        }
        // write out the cluster to the fragSize cluster file.
        // this should be it's own little function.
        if (currCluster.size() > 1) 
        {
            clusterId++;
            WritePosClusterToFile(clusters, currCluster, clusterId);
        } 
        
        // reset for the next cluster
        currCluster.clear();
        prevLine = currLine;
        currCluster.push_back(prevLine);
    }
}


void HydraPE::FindInterClusters(ifstream *sortedMappings, ofstream *clusters, int &clusterId) 
{
    PAIR prevLine;
    PAIR currLine;
    vector<PAIR> currCluster;
    
    *sortedMappings >> prevLine;
    currCluster.push_back(prevLine);
    while (*sortedMappings >> currLine) 
    {
        while (((currLine.core.start1 - prevLine.core.end1) <= 2000) &&
              sortedMappings->good())
        {
            currCluster.push_back(currLine);
            prevLine = currLine;
            *sortedMappings >> currLine;
        }
        // write out the cluster to the fragSize cluster file.
        // this should be it's own little function.
        if (currCluster.size() > 1) 
        {
            clusterId++;
            WritePosClusterToFile(clusters, currCluster, clusterId);
        } 
        
        // reset for the next cluster
        currCluster.clear();
        prevLine = currLine;
        currCluster.push_back(prevLine);
    }
}


/********************************************************************
 GOAL: Examine each size & pos cluster and choose the best "seed"
       mapping from which to assemble the breakpoint cluster.  
*********************************************************************/
void HydraPE::AssembleClusters() {

    cerr << "Assembling raw breakpoint clusters." << endl;            
        
    // assemble each cluster file that remained
    // after size and position sorting.
    int masterClusterId = 0;
    
    // sort the cluster files alphabetically for consistency with threads.
    sort(_posClusterFiles.begin(), _posClusterFiles.end());
    
    //#pragma omp parallel for
    for (int clusterNum = 0; clusterNum < _posClusterFiles.size(); ++clusterNum) {
        
        // open the current cluster file.
        string clusterFile = _posClusterFiles[clusterNum];
        ifstream cluster(clusterFile.c_str(), ios::in);
        
        string outfile = clusterFile + ".assembled";
        string puntfile = clusterFile + ".punted";
        ofstream *assembledClusters = new ofstream(outfile.c_str(), ios::out);
        ofstream *puntedClusters = new ofstream(puntfile.c_str(), ios::out);
        // load the mappings for this cluster into a vector in memory
        vector<PAIR> clusterMappings;
        PAIR mapping;
        int prevCluster = -1;
        while (cluster >> mapping) {
            if (mapping.clusterId != prevCluster && prevCluster > 0) {
                cout << "\tAssembling cluster " <<  prevCluster << " from: " << clusterFile << ".  N= " << clusterMappings.size() << endl;
                if (clusterMappings.size() <= 100000) {
                    Assemble(assembledClusters, clusterMappings, masterClusterId);
                }
                // too many mappings. save, but skip it.
                else {
                    cout << "\t\tSkipping " <<  prevCluster << " b/c it is too large. Locus and num mappings follows." << endl;
                    cout << "\t\t" << clusterMappings[0].core.chrom1 << "\t"
                                   << clusterMappings[0].core.start1 << "\t"
                                   << clusterMappings[clusterMappings.size()-1].core.start1 << "\t"
                                   << clusterMappings.size()
                                   << endl;
                    // write this cluster to the punted file.
                    for (int c = 0; c < clusterMappings.size(); c++) {
                        *puntedClusters << clusterMappings[c];
                    }
                }
                clusterMappings.clear();
                clusterMappings.push_back(mapping);
            }
            else {
                clusterMappings.push_back(mapping);
            }
            prevCluster = mapping.clusterId;
        }
        // process the last cluster
        Assemble(assembledClusters, clusterMappings, masterClusterId);
        // close the current cluster file.
        cluster.close();
        // close the assemble clusters file.
        assembledClusters->close();
        puntedClusters->close();
        cout << "FINISHED assembling clusters from " <<  clusterFile << "." << endl;
    }
}


void HydraPE::Assemble(ofstream *assembledClusters, vector<PAIR> &clusterMappings, int &clusterNum) {
    //***************************************************************
    // Phase 1.
    // Compute the support for each read in this cluster,
    // regardless of the mapping class.
    //***************************************************************
    int  nonOverlap;
    int  lengthDiff;
    int  support;
    long totalSupport;
    
    for (size_t i = 0; i < clusterMappings.size(); i++) {

        support = 0;
        for (size_t j = 0; j < clusterMappings.size(); j++) {

            if (j != i) {
                // changed to 2*
                int  aLengthDev         = GetLengthDeviation(clusterMappings[i]);
                int  bLengthDev         = GetLengthDeviation(clusterMappings[j]);
                int  aSpanDev           = GetSpanDeviation(clusterMappings[i]);
                int  bSpanDev           = GetSpanDeviation(clusterMappings[j]);
                
                bool haveLengthSupport  = doMultiLibraryLengthsSupportOneAnother(clusterMappings[i],
                                                                                 clusterMappings[j],
                                                                                 aLengthDev,
                                                                                 bLengthDev);
                bool haveSpanSupport    = doMultiLibrarySpansSupportOneAnother(clusterMappings[i],
                                                                               clusterMappings[j],
                                                                               aSpanDev,
                                                                               bSpanDev);
                if ( haveLengthSupport && haveSpanSupport ) 
                    support++;  
            }
        }
        // update the support for this mapping
        clusterMappings[i].support = support;
        totalSupport += support;
    }
            
    //***************************************************************
    // Phase 2.
    // - Sort the mappings by the support computed in Phase 1.
    // -  
    //*******************************************************************
    sort(clusterMappings.begin(), clusterMappings.end(), byTypeAndSupport);
    

    // tracks whether or not a mapping has already been used.
    std::map<string, bool> readAlreadyUsed;  
    int pass = 0;

    while(clusterMappings.size() > 1) {
        vector<PAIR> cluster;
        vector<int> cmIndexes;
        if (clusterMappings[0].support > 0) {

            // the list of pairs that support the seed...aka a contig
            cluster.push_back(clusterMappings[0]);
            cmIndexes.push_back(0);

            // remove the seed from the pool of pairs
            clusterMappings[0].used = true;
            readAlreadyUsed[clusterMappings[0].core.readId] = true;

            // compute the non-overlap between all the other reads and the seed
            for (int i = 0; i < clusterMappings.size(); i++) {
                clusterMappings[i].nonOverlap = getNonOverlap(clusterMappings[i], cluster[0]);
            }

            // sort the reads by non-overlap relative to the seed.
            sort(clusterMappings.begin(),clusterMappings.end(),byNonOverlap);
            
            // Now, loop through the pairs in order of how well they overlap
            // with the seed.  That is, in ascending order of "non-overlap".

            for (int currPair = 1; currPair < clusterMappings.size(); currPair++) {

                bool supportsAll = true;
                // Make sure the currPair is in support with ALL of the other pairs in the contig
                int clusterSize = cluster.size();
                
                bool haveLengthSupport, haveSpanSupport;
                for(int clusterPair = 0; clusterPair < clusterSize; clusterPair++) {
                    int  aLengthDev         = GetLengthDeviation(cluster[clusterPair]);
                    int  bLengthDev         = GetLengthDeviation(clusterMappings[currPair]);
                    int  aSpanDev           = GetSpanDeviation(cluster[clusterPair]);
                    int  bSpanDev           = GetSpanDeviation(clusterMappings[currPair]);
                    
                    haveLengthSupport = doMultiLibraryLengthsSupportOneAnother(cluster[clusterPair], 
                                                                                    clusterMappings[currPair], 
                                                                                    aLengthDev, 
                                                                                    bLengthDev);
                    haveSpanSupport   = doMultiLibrarySpansSupportOneAnother(cluster[clusterPair],
                                                                                  clusterMappings[currPair], 
                                                                                  aSpanDev,
                                                                                  bSpanDev);
                                                                                  
                    if ( (haveLengthSupport == false) || (haveSpanSupport == false)) {
                        supportsAll = false;
                        break;  // no need to continue testing if it doesn't support >=1 other mapping
                    }
                }
                // If the current pair is in support with ALL the other
                // pairs in the contig, then add it to the contig and remove it from
                // the working list of pairs.
                if (supportsAll == true) {
                    cluster.push_back(clusterMappings[currPair]);
                    cmIndexes.push_back(currPair);
                }
            }
    

            // add this contig to the list of all contigs
            // IF there is more than just the seed
            if (cluster.size() >= this->minSupport) {
                clusterNum++;
                for (int c = 0; c < cluster.size(); c++) {
                    // remove each read from the "pool" of mappings 
                    // (cmIndexes is used to get the original index for this mappings)
                    clusterMappings[cmIndexes[c]].used      = true;
                    cluster[c].clusterId = clusterNum; // convert to 1-based
                    *assembledClusters << cluster[c];
                    totalSupport           -= clusterMappings[c].support;
                }
                cluster.clear();
                cmIndexes.clear();
            }
            else {
                cluster.clear();
                cmIndexes.clear();
            }
            // sort the reads so the unused mappings are at the "top"
            sort(clusterMappings.begin(),clusterMappings.end(),byUsed);

            // delete the used mappings
            vector<PAIR>::iterator it = clusterMappings.begin();
            while ((it != clusterMappings.end()) && !(it->used)) {
                it++;
            }
            clusterMappings.erase(it, clusterMappings.end());

            // resort the reads by type and support for the next seed search.
            sort(clusterMappings.begin(),clusterMappings.end(),byTypeAndSupport);       
        }
        else {
            clusterMappings.erase(clusterMappings.begin());
        }   
        pass++;
    }
}


//***************************************************
// Method: SummarizeContigs
// Purpose:
//
// Outline:
//***************************************************
void HydraPE::ReportSVCalls(string &outStub) {

    string allFile    = outStub;
    string finalFile  = outStub;
    string detailFile = outStub;
    
    detailFile.append(".detail");
    finalFile.append(".final");
    allFile.append(".all");
    
    // open the all contig file for writing
    ofstream all(allFile.c_str(), ios::out);
    if ( !all ) {
        cerr << "Error: The file of _all_ breakpoints (" << allFile << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    
    // open the finall contig file for writing
    ofstream final(finalFile.c_str(), ios::out);
    if ( !final ) {
        cerr << "Error: The file of _final_ breakpoints (" << finalFile << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }
    
    // open the detail contig file for writing
    ofstream detail(detailFile.c_str(), ios::out);
    if ( !detail ) {
        cerr << "Error: The file of breakpoints details (" << detailFile << ") could not be opened.  Exiting!" << endl;
        exit (1);
    }

    int contigId = 0;                       // incremental id for each event that has the required support.
    int numMappings;                        // number of mappings in contig
    
    int minStart1, minStart2;               // used to track the "footprint" of contig
    int maxEnd1, maxEnd2;
    
    int totalMM1, totalMM2;                 // summary info
    int numUniqueMappers, numAnchoredMappers, numMultipleMappers;
    int totalMappings1, totalMappings2;
    int totalQual1, totalQual2;
    int meanMappings1, meanMappings2;
    int meanQual1, meanQual2;
    float meanMM1, meanMM2;
            
    double allWeightedSupport;
    int finalSupport;
    double finalWeightedSupport;
            
    // loop through all the contigs that were found 
    for (int j = 0; j < this->contigs.size(); j++) {

        vector<PAIR> contigMappings  = contigs[j];              // get the mappings for this contig
        numMappings                  = contigMappings.size();   // get the number of mappings for this contig
        
        // The "footprint" of the breakpoint
        minStart1 = minStart2 = INT_MAX;    
        maxEnd1   = maxEnd2   = 0;
        
        totalMM1         = totalMM2 = 0;
        numUniqueMappers = numAnchoredMappers = numMultipleMappers = 0;
        totalMappings1   = totalMappings2 = 0;
        totalQual1       = totalQual2 = 0;
        meanMappings1    = meanMappings2 = 0;
        meanQual1        = meanQual2 = 0;
        meanMM1          = meanMM2 = 0.0;
        
        allWeightedSupport   = 0.0;
        finalSupport         = 0;
        finalWeightedSupport = 0.0;
        
        int numAllUniques    = 0;
        int numFinalUniques  = 0;
        int contigSize;  // what is the "span" of the breakpoint
        
        // get the contig size and "footprints" based on the min and max mapping coordinates
        // if there is final support, only use the mappings that were _included_ in the cluster
        // otherwise, define the footprint based on all the mappings that were once in the cluster
        if (hasFinalSupport(contigMappings) == true) {
            contigSize = getClusterSizeFinal(contigMappings, minStart1, minStart2, maxEnd1, maxEnd2);
        }
        else {
            contigSize = getClusterSizeAll(contigMappings, minStart1, minStart2, maxEnd1, maxEnd2);         
        } 

        // (1)  Get _all_ of the unique ids in this cluster
        getNumUniquePairs(contigMappings, numFinalUniques, numAllUniques);
        
        // (2)  Get the total edit distance on each end of this cluster
        getTotalEditDistance(contigMappings, totalMM1, totalMM2);
        
        // (3)  Get the total number of mappings on each end of this cluster
        getTotalNumMappings(contigMappings, totalMappings1, totalMappings2);
        
        // (4)  Tally the final, finalWeighted and allWeightedSupport for this cluster
        computeSupport(contigMappings, finalSupport, finalWeightedSupport, allWeightedSupport, 
                        numUniqueMappers, numAnchoredMappers, numMultipleMappers);

        // calculate quality statistics for each contig.
        if (finalSupport) {
            meanMappings1 = totalMappings1 / finalSupport;
            meanMappings2 = totalMappings2 / finalSupport;
            meanQual1     = totalQual1 / finalSupport;
            meanQual2     = totalQual2 / finalSupport; 
            meanMM1       = float(totalMM1) / float(finalSupport);
            meanMM2       = float(totalMM2) / float(finalSupport);
        }

        // report this contigs as longs as there is at least one uniq
        // pair in the contig.
        if (allWeightedSupport > 0) {
            
            contigId++;
                        
            // BEDPE-style
            all << contigMappings[0].core.chrom1 << "\t" << minStart1 << "\t" << maxEnd1 << "\t" <<
                   contigMappings[0].core.chrom2 << "\t" << minStart2 << "\t" << maxEnd2 << "\t" <<
                   contigId << "\t" << numFinalUniques << "\t" << contigMappings[0].core.strand1 << "\t" << 
                   contigMappings[0].core.strand2 << "\t" << meanMM1 << "\t" << meanMM2 << "\t" << 
                   meanMappings1 << "\t" << meanMappings2 << "\t" << contigSize << "\t" << numMappings << "\t" << 
                   allWeightedSupport << "\t" << finalSupport << "\t" << finalWeightedSupport << "\t" << 
                   numUniqueMappers << "\t" << numAnchoredMappers << "\t" << numMultipleMappers << "\t" << endl;
                    
            if (numFinalUniques >= this->minSupport) {
                // BEDPE-style
                final << contigMappings[0].core.chrom1 << "\t" << minStart1 << "\t" << maxEnd1 << "\t" <<
                       contigMappings[0].core.chrom2 << "\t" << minStart2 << "\t" << maxEnd2 << "\t" <<
                       contigId << "\t" << numFinalUniques << "\t" << contigMappings[0].core.strand1 << "\t" << 
                       contigMappings[0].core.strand2 << "\t" << meanMM1 << "\t" << meanMM2 << "\t" << 
                       meanMappings1 << "\t" << meanMappings2 << "\t" << contigSize << "\t" << numMappings << "\t" << 
                       allWeightedSupport << "\t" << finalSupport << "\t" << finalWeightedSupport << "\t" << 
                       numUniqueMappers << "\t" << numAnchoredMappers << "\t" << numMultipleMappers << "\t" << endl;
            }   

            // write the read/mapping information for this contig to the contig detail file.
            for (int i = 0; i < numMappings; i++) {
                if (contigMappings[i].include) {
                    detail << contigMappings[i].core.chrom1 << "\t" << contigMappings[i].core.start1 << "\t" << contigMappings[i].core.end1 << "\t"
                           << contigMappings[i].core.chrom2 << "\t" << contigMappings[i].core.start2 << "\t" << contigMappings[i].core.end2 << "\t"
                           << contigMappings[i].core.readId << "\t" << contigMappings[i].core.mate1 << "\t" << contigMappings[i].core.strand1 << "\t"
                           << contigMappings[i].core.strand2 << "\t" << contigMappings[i].core.edit1 << "\t" << contigMappings[i].core.edit2 << "\t"
                           << contigMappings[i].core.mappings1 << "\t" << contigMappings[i].core.mappings2 << "\t" 
                           << int(contigMappings[i].mappingType) << "\t" << "Y" << "\t" << contigId << endl;                                        
                }
                else {
                    detail << contigMappings[i].core.chrom1 << "\t" << contigMappings[i].core.start1 << "\t" << contigMappings[i].core.end1 << "\t"
                           << contigMappings[i].core.chrom2 << "\t" << contigMappings[i].core.start2 << "\t" << contigMappings[i].core.end2 << "\t"
                           << contigMappings[i].core.readId << "\t" << contigMappings[i].core.mate1 << "\t" << contigMappings[i].core.strand1 << "\t"
                           << contigMappings[i].core.strand2 << "\t" << contigMappings[i].core.edit1 << "\t" << contigMappings[i].core.edit2 << "\t"
                           << contigMappings[i].core.mappings1 << "\t" << contigMappings[i].core.mappings2 << "\t" 
                           << int(contigMappings[i].mappingType) << "\t" << "N" << "\t" << contigId << endl;
                }
            }
        }
    }
    all.close();
    final.close();
    detail.close();
}



//***********************************************************
// Method: CullMappingsByMisMatches
// Purpose: To reduce the mappings for a given pair
//          to those with the least mismatches.
//
// Outline: 
//  1.  Find the mapping(s) with the least mismatches.
//  2.  Sort the mappings so the least MM set is on "top".   
//  3.  Ditch all but the least MM set.
//***********************************************************
void HydraPE::CullMappingsByMisMatches(pairVector &pairMappings) {

    // set of unique mappings for each mate in a pair
    std::map<BED3, bool> end1Mappings, end2Mappings;
        
    unsigned short minMM = USHRT_MAX;
    int totalMM;
    
    // loop through the mappings and track the unique 
    // mappings on each end as well as the minimum edit
    // distance observed for any one mapping.  this minimum
    // is the standard for all other mappings.  
    pairVector::const_iterator mapIter = pairMappings.begin();
    pairVector::const_iterator mapEnd = pairMappings.end();
    for (; mapIter != mapEnd; ++mapIter) {      
        if (mapIter->core.mate1 == 1) {
            BED3 end1Mapping(mapIter->core.chrom1, mapIter->core.start1, mapIter->core.end1);
            BED3 end2Mapping(mapIter->core.chrom2, mapIter->core.start2, mapIter->core.end2);
            end1Mappings[end1Mapping] = true;
            end2Mappings[end2Mapping] = true;           
        }
        else {
            BED3 end2Mapping(mapIter->core.chrom1, mapIter->core.start1, mapIter->core.end1);
            BED3 end1Mapping(mapIter->core.chrom2, mapIter->core.start2, mapIter->core.end2);
            end1Mappings[end1Mapping] = true;
            end2Mappings[end2Mapping] = true;       
        }
        
        // does this mapping have the least mismatches?
        totalMM = getTotalMM(*mapIter);
        if (totalMM < minMM)
            minMM = totalMM;
    }   

    // Figure out what type of pair this is based on the alignments.
    short mapType;
    int end1Size = end1Mappings.size();
    int end2Size = end2Mappings.size();
    
    if ((end1Size == 1) && (end2Size == 1))
        mapType = UNIQ_TYPE;
    else if ((end1Size == 1) || (end2Size == 1))
        mapType = ANCH_TYPE;
    else
        mapType = MULT_TYPE;

    // sort the reads so the mappings with the least mismatches are at the "top"                                                                                                                                                    
    sort(pairMappings.begin(), pairMappings.end(), byTotalMM);
                                                                                                                                        
    // default to assuming we want to use just the best mappings
    // "withinBest" = allow up to editBeyondBest edits worse than minMM
    // "all" = we don't want to delete anything.
    int editDistanceCutoff = minMM;
    if (this->mappingUsage == "withinBest") 
        editDistanceCutoff = minMM + editBeyondBest;
    else if (this->mappingUsage == "all")
        editDistanceCutoff = INT_MAX;
    
    vector<PAIR>::iterator mappingsIter = pairMappings.begin();
    while ((mappingsIter != pairMappings.end()) && (getTotalMM(*mappingsIter) <= editDistanceCutoff)) {
        mappingsIter->mappingType = mapType; // set the mapping type for all of the mappings within the edit distance cutoff
        mappingsIter++;
    }

    // ditch all but the least mm set.
    pairMappings.erase(mappingsIter, pairMappings.end());
    
    // Now, store the number of mappings observed on each end of this pair
    pairVector::iterator pairIter = pairMappings.begin();
    pairVector::iterator pairEnd  = pairMappings.end();
    for (; pairIter != pairEnd; ++pairIter) {
        pairIter->core.mappings1 = end1Size;
        pairIter->core.mappings2 = end2Size;
    }       
}
