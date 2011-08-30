/********************************************************************************
Program: 		Sort.cpp
Description: 	Ancillary comparison functions for sorting Hydra data structures.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*********************************************************************************/
#include "Sort.h"


/******************************************************************************
   Sorting comparison functions.
*******************************************************************************/
bool bySupport(const PAIR &a, const PAIR &b) {
	return (a.support < b.support);
}


bool byNonOverlap(const PAIR &a, const PAIR &b) {
	return (a.nonOverlap < b.nonOverlap);
}


bool byMappingType(const PAIR &a, const PAIR &b) {
	return (a.mappingType < b.mappingType);
}


bool byTypeAndSupport(const PAIR &a, const PAIR &b) {

	// first, sort by Mapping Type (asc)
	if      (a.mappingType < b.mappingType) return true;
	else if (a.mappingType > b.mappingType) return false;

	// second, within a mapping type, sort by support (desc)
	if (a.support > b.support) return true;
	else return false;
}


bool byUsed(const PAIR &a, const PAIR &b) {
	return (a.used < b.used);
}


bool byRead(const PAIR &a, const PAIR &b) {
	return (a.core.readId < b.core.readId);
}


bool byReadAndType(const PAIR &a, const PAIR &b) {

	// first, sort by read (asc)
	if      (a.core.readId < b.core.readId) return true;
	else if (a.core.readId > b.core.readId) return false;

	// second, within a read, sort by mapping type (asc)
	if (a.mappingType < b.mappingType) return true;
	else return false;
}


bool byTotalMM(const PAIR &a, const PAIR &b) {
	return (getTotalMM(a) < getTotalMM(b));
}


bool byWeightedSupportAndTotalMM(const IN_CLUSTER &a, const IN_CLUSTER &b) {

	// first, sort by weighted support (desc)
	if      (a.weightedSupport > b.weightedSupport) return true;
	else if (a.weightedSupport < b.weightedSupport) return false;

	// second, sort by totalMM (asc)
	if (a.totalMM < b.totalMM) return true;
	else return false;
}


bool byWeightedSupportAndSize(const IN_CLUSTER &a, const IN_CLUSTER &b) {

	// first, sort by weighted support (desc)
	if      (a.weightedSupport > b.weightedSupport) return true;
	else if (a.weightedSupport < b.weightedSupport) return false;

	// second, sort by size (asc)
	if (a.contigSize < b.contigSize) return true;
	else return false;
}


bool byFragmentSize(const PAIR &a, const PAIR  &b) {
	return (getFragSize(a.core) < getFragSize(b.core));
}

// New for multiple library support.
bool byDiffFromExpectedSize(const PAIR &a, const PAIR  &b) {
	return (a.diffFromExpectedSize < b.diffFromExpectedSize);
}


bool byStart1Asc(const PAIR &a, const PAIR &b) {
	return (a.core.start1 < b.core.start1);
}


bool byStart1Desc(const PAIR &a, const PAIR &b) {
	return (a.core.start1 >= b.core.start1);
}