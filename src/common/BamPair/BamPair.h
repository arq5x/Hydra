/****************************************************************************
   BamPair.h (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.

   Object for storing the BamAlignments (N=2) for 
   a readpair mapping combination
****************************************************************************/

#include "BamAux.h"
using namespace BamTools;

#include <iostream>
#include <stdlib.h>
using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BamPair {

public:

	// constructor 
	BamPair(const BamAlignment &end1, const BamAlignment &end2);
	// destructor
	~BamPair(void);
	// copy constructor
	BamPair(const BamPair &b);
	// returns true if valid ends of same pair.
	bool SetPairEnds();
	// return whether or not the pair is intra-chromosomal
	bool IsIntraChromosomal () const;
	// return the size (i.e. mapping distance) of the pair
	int GetSize () const;
	
	void SetPairAsConcordant();
    void SetPairAsDiscordant();

    BamAlignment GetFirstEndOfPair() const;
    BamAlignment GetSecondEndOfPair() const;
	
	// print the pair in Hydra format.
	void PrintPairForHydra(const RefVector &refs, int numEnd1Mappings, int numEnd2Mappings);
		
private:
	BamAlignment end1;
	BamAlignment end2;	
	

};
