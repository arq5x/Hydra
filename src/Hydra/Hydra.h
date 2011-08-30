#ifndef HYDRA_H
#define HYDRA_H

#include <algorithm>
#include <numeric>    // for the accumulate algorithm  
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <map>
#include <vector>
#include <list>

#include "DnaLibrary.h"


using namespace std;


// prevent duplicate includes of this header.
#pragma once


// Constants
const double UNIQ_WEIGHT = 1.0;
const double ANCH_WEIGHT = 0.5;
const double MULT_WEIGHT = 0.1;

const short UNIQ_TYPE = 1;
const short ANCH_TYPE = 2;
const short MULT_TYPE = 3;

const int BUFFER_SIZE = 5000000;

// Main stuctur for reading input BEDPE records.
struct CORE_PAIR {
    string readId;                		// the read id.

	string chrom1, chrom2;     			// the mapped chrom for read1, read2
	short mate1, mate2;                 // which end of the mate-pair
	int start1, start2;           		// the mapped start position for read1, read2
	int end1, end2;           			// the mapped end position for read1, read2
	string strand1, strand2;   			// the mapping orientation for read1, read2
    int whichMateIsBlock1;              // 1 if the first coordinate block is from end 1 of the pair. 2 otherwise.
	short edit1, edit2;                 // the number of mismatches for read1, read2
	int mappings1, mappings2;           // the number of mappings observed on each end.
	
	friend istream& operator>>(istream &is, CORE_PAIR &p) {
        is  >> p.chrom1    >> p.start1            >> p.end1
			>> p.chrom2    >> p.start2            >> p.end2
			>> p.readId    >> p.whichMateIsBlock1
			>> p.strand1   >> p.strand2
            >> p.edit1     >> p.edit2;
        return is;
    }
};


// Main structure to store a specific mapping for a paired-end alignment
struct PAIR {

    CORE_PAIR core;
    
	short mappingType;        		    // the mapping type: 1 = dis_unique, 2 = anchor, 3 = multiple, 4 = orphan
	unsigned int nonOverlap;            // this is a temporary variable used to compute the 
									    // number of bases in this pair the DO NOT overlap
								        // with the current contig seed that is being evaluated.
	unsigned int support;               // the number of reads in support of this read
    int clusterId;
	bool used;
	bool include;
    unsigned short fileNum;             // the index of the discordant filenum from the command line (0-based)
    
    unsigned int diffFromExpectedSize;
    
    friend ostream& operator<<(ostream &os, const PAIR &p) {
        os  << p.core.chrom1          << "\t" << p.core.start1            << "\t" << p.core.end1 << "\t"
	        << p.core.chrom2          << "\t" << p.core.start2            << "\t" << p.core.end2 << "\t"
			<< p.core.readId          << "\t" << p.core.whichMateIsBlock1 << "\t"
			<< p.core.strand1         << "\t" << p.core.strand2           << "\t"
            << p.core.edit1           << "\t" << p.core.edit2             << "\t"
            << p.core.mappings1       << "\t" << p.core.mappings2         << "\t"
            << p.mappingType          << "\t"
            << p.diffFromExpectedSize << "\t"
            << p.support              << "\t"
            << p.clusterId            << "\t"
            << p.used                 << "\t"
            << p.include              << "\t"
            << p.fileNum              << "\t"
            << endl;
        return os;
    }
    
    friend istream& operator>>(istream &is, PAIR &p) {
        is  >> p.core.chrom1 >> p.core.start1  >> p.core.end1
			>> p.core.chrom2 >> p.core.start2  >> p.core.end2
			>> p.core.readId >> p.core.whichMateIsBlock1
			>> p.core.strand1 >> p.core.strand2
            >> p.core.edit1   >> p.core.edit2
            >> p.core.mappings1 >> p.core.mappings2
            >> p.mappingType
            >> p.diffFromExpectedSize
            >> p.support
            >> p.clusterId
            >> p.used
            >> p.include
            >> p.fileNum;
        return is;
    }
};


// structure for the coordinates of an alignment
// Using this struct as a key to a mathis->  Thus, need the "<" overload so that
// the map can compare BED3.
struct BED3 {
	
	// data
	string chrom;
	int start;
	int end;
	
	// constructor
	BED3 (const string c, const int s, const int e) :
		chrom(c), 
		start(s), 
		end(e) 
	{}
	
	// comparison operator for maps keyed on this structure
	bool operator < (const BED3 &a) const
	{ 
		// first, sort by chrom (asc)
		if (chrom < a.chrom)      return true;
		else if (chrom > a.chrom) return false;

		// second, within a chrom sort by start
		if (start < a.start)      return true;
		else if (start > a.start) return false;
		
		if (end < a.end)      return true;
		else return false;
	}
};


// structure to track the "contigs" that a given mapping 
// has been assigned to. 
struct IN_CLUSTER {
	int contig;
	bool include;
	double weightedSupport;
	int totalMM;
	int contigSize;
	bool unlinked;
};


// Using this struct as a key to a mathis->  Thus, need the "<" overload so that
// the map can compare CHROMS_AND_STRANDS.
struct CHROMS_AND_STRANDS {
	// data
	string chrom1;
	string chrom2;
	string strand1;
	string strand2;
	
	// constructor
	CHROMS_AND_STRANDS (const string chrom1, const string chrom2,
						const string strand1, const string strand2) :
		chrom1(chrom1),
		chrom2(chrom2),
		strand1(strand1),
		strand2(strand2) 
	{}
	
	// comparison operator for maps keyed on this structure
	bool operator < (const CHROMS_AND_STRANDS &a) const
	{ 
		// first, sort by chrom1 (asc)
		if (chrom1 < a.chrom1)      return true;
		else if (chrom1 > a.chrom1) return false;

		// second, sort by chrom2 (asc)
		if (chrom2 < a.chrom2)      return true;
		else if (chrom2 > a.chrom2) return false;
		
		// third, sort by strand1 (asc)
		if (strand1 < a.strand1)      return true;
		else if (strand1 > a.strand1) return false;
		
		// lastly, sort by strand2 (asc)
		if (strand2 < a.strand2) return true;
		else return false;
	}
};



// TYPEDEFS for common data structures


typedef map<CHROMS_AND_STRANDS, vector< vector<PAIR> >, std::less<CHROMS_AND_STRANDS> > putativeContigMap;

typedef multimap<int, vector<PAIR>, std::less<int> > circosContigMap;

typedef map<string, vector<IN_CLUSTER>, std::less<string> > read2ContigsMap;

typedef map<string, vector<PAIR>, std::less<string> > read2MappingsMap; 

typedef map<string, ostream*>::const_iterator outFileNameToPointerMap;

// vector of paired-end mappings
typedef vector<PAIR> pairVector;

// map of paired-end mappings key by chrom-chrom-strand-strand
typedef map<CHROMS_AND_STRANDS, vector<PAIR>, std::less<CHROMS_AND_STRANDS> > pairMap;

// vector of contigs containing all mappings in each contig
typedef vector< vector<PAIR> > contigVector;


//*************************************************
// Template Functions
//*************************************************

// templated function to convert objects to strings
template <typename T>
std::string ToString(const T & value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}



class HydraPE {

public:

	//************************************************
	// Public methods and elements
	//************************************************	
	
	// constructor 
	HydraPE(vector<DNALIB> libraries, int minSupport, 
	        int maxLinkedDistance, bool ignoreSize, 
	        bool lumpInversions, string mappingUsage, int editBeyondBest, int memory);

	// destructor
	~HydraPE(void);

    // new stuff.
	void Tokenize(const string &, vector<string> &, const string &);  	// tokenize a line into a vector of strings
    void SetExpectedSizes(const vector<int> &expectedSizes);
    void SetVariances(const vector<int> &variances);
    void SetDegreesOfVariance(const vector<int> &degreesOfVariance);

    // mapping routing
	void RouteDiscordantMappings();				                        // route mappings to files based on chr/chr/strand/strand
    void LoadRoutedFile(const string &routedFile);                      // start with a single already routed file.
    void LoadRoutedFileList(const string &routedFiles);                 // start with list of already routed files.
    void LoadPosSortedFileList(const string &sortedFiles);
    void LoadPosClusterFileList(const string &posClusterFiles);         // start with already created posCluster files.

    // sort mapping files
	void SortFragments();												// find clusters of discordant mappings
    // void SortAllMasterFilesByFragSize();
    //void SortAllMasterFilesByPosition();
    void SortAllMasterFilesByPosition_New();

    // find clusters in sorted files
	// void FindFragSizeClusters();
    // void FindPositionClusters();
    void FindPositionClusters_New();
    
    // assemble frag/pos clusters into putative breakpoint calls.
    void AssembleClusters(int maxMappings);
    
    //void AllowOneClusterPerPair();
    
    // methods for removing unnecessary files.
    void RemoveMasterChromStrandFiles();
    void RemoveFragSizeSortedFiles();
    void RemovePositionSortedFiles();
    void RemoveFragSizeClusterFiles();
    void RemovePositionClusterFiles();
    
    // old stuff.
	void AssembleContigs();												// assemble clusters of discordant mappings
	void AllowOneContigPerRead();										// find the most likely contig for each pair
	void ReportSVCalls(string &);										// report the final SV calls to output files.


private:

	//************************************************
	// Private methods and elements
	//************************************************
    vector<string> _discordantFiles;                  // the original discordant files from CL
    vector<int>    _expectedSizes;
    vector<int>    _variances;
    vector<int>    _degreesOfVariance;
    int            _memory;
    int            _maxMappings;
    
    vector<DNALIB> _libraries;

	// the intermediate files
    map<string, ostream*> _masterChromStrandFiles;     // the name and ostream of the chr/chr/str/str files
    map<string, bool>     _masterChromStrandFileTypes; // is the file intra- (True) or inter- (False) chromosomal

	//vector<string>        _fragSortFiles;             // the names of the master files that are sorted by frag size
	//vector<string>        _fragClusterFiles;          // the name and ostream of frag size cluster files

	vector<string>        _posSortFiles;              // the names of the master files that are sorted by position
    vector<bool>          _posSortFileTypes;          // are the files intra- (True) or inter (False)
	vector<string>        _posClusterFiles;           // the name and ostream of position cluster files
    
	int lengthDev;
	int spanDev;
	int minSupport;
	int maxLinkedDistance;
	bool ignoreSize;
	bool lumpInversions;
	string mappingUsage;
	int editBeyondBest;
	
	// private data structures for main processing.
	pairMap mappingsByChromAndStrand;
	putativeContigMap putativeContigs;
	contigVector contigs;
	read2ContigsMap  _read2Clusters;
	read2MappingsMap read2Mappings;
	
	// private methods
	void CullMappingsByMisMatches(pairVector &mappings);	
	void CorrectMateOrder(CORE_PAIR &pair);
	void SwapEnds(CORE_PAIR &pair);
	void AddMappingsToMasterMap(const pairVector &mappings);
	ostream* GetMasterFileStream(const string &masterFileName, bool isIntra);
    bool IsInversionMapping(const CORE_PAIR &pair);
	void AddMappingsToMasterFile(const vector<pairVector> &readMappings);
    int  GetLengthDeviation(const PAIR &pair);
    int  GetSpanDeviation(const PAIR &pair);
    void WriteFragClusterToFile(ofstream *out, vector<PAIR> &cluster, int clusterId);
    void WritePosClusterToFile(ofstream *out, vector<PAIR> &cluster, int clusterId);
    void Assemble(ofstream *assembledClusters, vector<PAIR> &clusterMappings, int &clusterNum);
    void ExamineFragmentCluster(vector<PAIR> &fragClusterMappings, int &clusterId, ofstream *output);
    void RouteFile(const string &file, int fileNum);
    void FindIntraClusters(ifstream *sortedMappings, ofstream *clusters, int &clusterId);
    void FindInterClusters(ifstream *sortedMappings, ofstream *clusters, int &clusterId);    
};

#endif
