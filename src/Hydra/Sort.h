/********************************************************************************
Program: 		Sort.h
Description: 	Ancillary comparison functions for sorting Hydra data structures.
					
Author:    		Aaron Quinlan, Ph.D
				University of Virginia
				aaronquinlan@gmail.com
*********************************************************************************/
#include "Ancillary.h"

/******************************************************************************
   Sorting comparison functions.
*******************************************************************************/
bool bySupport (const PAIR &a, const PAIR &b);
bool byNonOverlap (const PAIR &a, const PAIR &b);
bool byMappingType (const PAIR &a, const PAIR &b);
bool byTypeAndSupport (const PAIR &a, const PAIR &b);
bool byUsed (const PAIR &a, const PAIR &b);
bool byRead (const PAIR &a, const PAIR &b);
bool byReadAndType (const PAIR &a, const PAIR &b);
bool byTotalMM (const PAIR &a, const PAIR &b);
bool byWeightedSupportAndTotalMM (const IN_CLUSTER &a, const IN_CLUSTER &b);
bool byWeightedSupportAndSize (const IN_CLUSTER &a, const IN_CLUSTER &b);
bool byFragmentSize (const PAIR &a, const PAIR  &b);
bool byDiffFromExpectedSize(const PAIR &a, const PAIR  &b);
bool byStart1Asc (const PAIR &a, const PAIR &b);
bool byStart1Desc (const PAIR &a, const PAIR &b);