/********************************************************************************
Program: 		Ancillary.h
Description: 	Ancillary "helper" functions for Hydra data structures.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*********************************************************************************/
#include "Hydra.h"

// Return the "inner span" of a mapping
int    getFragSize (const CORE_PAIR &pair);

bool   doSpansSupportCommonBreakpoint (const PAIR &a, const PAIR &b);

int    getNonOverlap (const PAIR &a, const PAIR &b);

bool   doLengthsSupportOneAnother (const PAIR &a, const PAIR &b, int lengthDev);

bool   doMultiLibraryLengthsSupportOneAnother (const PAIR &a, const PAIR &b, int aLengthDev, int bLengthDev); 

bool   doSpansSupportOneAnother (const PAIR &a, const PAIR &b, int spanDev);

bool   doMultiLibrarySpansSupportOneAnother(const PAIR &a, const PAIR &b, int aSpanDev, int bSpanDev);

int    getTotalMM (const PAIR &pair);

int    getTotalMMAmongAllMappings(const pairVector &mappings);

double getTotalWeightedSupportAmongAllMappings(const pairVector &mappings);

// Return true if a cluster has at least one mapping that is "included"
bool   hasFinalSupport(const pairVector &mappings);

int    getClusterSizeFinal(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2);


int    getClusterSizeAll(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2);


bool   isVariantUnlinked(const pairVector &mappings, int maxDistance);


bool   areAnyContigsUnlinked(const vector<IN_CLUSTER> &contigs);


void   getNumUniquePairs(const pairVector &mappings, int &numFinalUniques, int &numAllUniques);


void   getTotalEditDistance(const pairVector &mappings, int &totalMM1, int &totalMM2);


void   getTotalNumMappings(const pairVector &mappings, int &totalMappings1, int &totalMappings2);
	

void   computeSupport(const pairVector &mappings, int &finalSupport, double &finalWeightedSupport, double &allWeightedSupport, 
                      int &numUniqueMappers, int &numAnchoredMappers, int &numMultipleMappers);

