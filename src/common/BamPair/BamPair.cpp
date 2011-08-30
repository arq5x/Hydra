/****************************************************************************
   BamPair.cpp (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.

   Object for storing the BamAlignments (N=2) for 
   a readpair mapping combination
****************************************************************************/
#include "BamPair.h"

using namespace std;


// constructor
BamPair::BamPair (const BamAlignment &end1, const BamAlignment &end2) 
	: end1(end1)
	, end2(end2)
{}


// destructor
BamPair::~BamPair(void) {}

// copy constructor
BamPair::BamPair(const BamPair &b) {
	end1 = b.end1;
	end2 = b.end2;
}

// set the alignments for each end of this BAM pair. 
bool BamPair::SetPairEnds () {
			
	if (end1.RefID == end2.RefID) {	
		if (end2.Position < end1.Position) {
			BamAlignment temp = end1;
			end1 = end2;
			end2 = temp;
            return true;
		}
        return false;
	}
	else {	
		if (end2.RefID < end1.RefID) {
			BamAlignment temp = end1;
			end1 = end2;
			end2 = temp;
            return true;
		}
        return false;
	}
}


bool BamPair::IsIntraChromosomal () const {
	return end1.RefID == end2.RefID;	
}


int BamPair::GetSize () const {
	return abs(end1.Position - end2.Position);	
}

void BamPair::SetPairAsConcordant() {
    end1.SetIsProperPair(true);
    end2.SetIsProperPair(true);
}

void BamPair::SetPairAsDiscordant() {
    end1.SetIsProperPair(false);
    end2.SetIsProperPair(false);    
}

BamAlignment BamPair::GetFirstEndOfPair() const {
    return end1;
}

BamAlignment BamPair::GetSecondEndOfPair() const {
    return end2;
}
    
/*
	Report the alignment in Hydra format

	Hydra Layout:
	1. read
	2. chrom1
	3. pos1
	4. strand1
	5. mm1
	6. chrom2
	7. pos2
	8. strand2
	9. mm2
*/ 
void BamPair::PrintPairForHydra(const RefVector &refs, int numEnd1Mappings, int numEnd2Mappings) {

	// resolve the strands for each end
	string strand1 = "+"; 
	string strand2 = "+";
	if (end1.IsReverseStrand()) strand1 = "-"; 
	if (end2.IsReverseStrand()) strand2 = "-";
	
	// resolve whic end of the pair the first block comes from
	unsigned short mate = 1;
    if (end1.IsFirstMate() == false) mate = 2;
		
	uint32_t ed1, ed2;
    
	// both ends are mapped
	if (end1.GetTag("NM", ed1) && end2.GetTag("NM", ed2)) {
		cout << refs.at(end1.RefID).RefName << "\t" << end1.Position << "\t" << end1.GetEndPosition(false) << "\t" 
		     << refs.at(end2.RefID).RefName << "\t" << end2.Position << "\t" << end2.GetEndPosition(false) << "\t"
		     << end1.Name << "\t" << mate << "\t" 
		     << strand1 << "\t" << strand2 << "\t"
		     << ed1 << "\t" << ed2 << "\t"
             << numEnd1Mappings << "\t" << numEnd2Mappings
             << endl;
	}
	// it's an orphan
	else if (end1.GetTag("NM", ed1)) {
		cout << refs.at(end1.RefID).RefName << "\t" << end1.Position << "\t" << end1.GetEndPosition(false) << "\t" 
		     << "." << "\t" << "-1" << "\t" << "-1" << "\t"
		     << end1.Name << "\t" << mate << "\t" 
		     << strand1 << "\t" << "." << "\t"
		     << ed1 << "\t" << "-1" << "\t"
             << numEnd1Mappings << "\t" << "-1"
             << endl;	    
	}
	// it's an orphan
	else if (end2.GetTag("NM", ed2)) {
		cout << refs.at(end2.RefID).RefName << "\t" << end2.Position << "\t" << end2.GetEndPosition(false) << "\t" 
		     << "." << "\t" << "-1" << "\t" << "-1" << "\t"
		     << end2.Name << "\t" << mate << "\t" 
		     << strand2 << "\t" << "." << "\t"
		     << ed2 << "\t" << "-1" << "\t"
             << numEnd2Mappings << "\t" << "-1"
             << endl;	    
	}
}

