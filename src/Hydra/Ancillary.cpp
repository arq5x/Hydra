/********************************************************************************
Program:        Ancillary.cpp
Description:    Ancillary helper functions for Hydra data structures.

Author:         Aaron Quinlan, Ph.D
                University of Virginia
                aaronquinlan@gmail.com
*********************************************************************************/
#include "Ancillary.h"


/******************************************************************************
   Utility functions.
*******************************************************************************/
int getFragSize (const CORE_PAIR &pair) {
    /* old, Contigger */
    //return (a.start2 - a.start1 + 1);
    return (pair.start2 - pair.end1 + 1);
}


bool doSpansSupportCommonBreakpoint (const PAIR &a, const PAIR &b) {
    /***************************************************************
     We want to enforce that two mappings have acceptable non-overlap
     AND that they support the same putative breakpoint region.

    For example, we want this:

         s1   e1                     s2    e2
    (A)   ======......................======
    (B)                   ======...............................=====
                          s1   e1                             s2   e2

    But NOT this:
         s1   e1                     s2    e2
    (A)   ======......................======
    (B)                          ======...............................=====
                                 s1   e1                             s2   e2
    ****************************************************************/

    if ((a.core.chrom1 == a.core.chrom2) && (b.core.chrom1 == b.core.chrom2)) {
        // int overlap = (min(a.core.end2, b.core.end2) - max(a.core.start1, b.core.start1));  // required outer-span overlap (wrong)
        int overlap = (min(a.core.start2, b.core.start2) - max(a.core.end1, b.core.end1));     // require inner-span overlap (correct)
        // if ( (overlap > 0) &&
        //      (a.core.end1 < b.core.start2) &&
        //      (b.core.end1 < a.core.start2))
        if (overlap > 0)
            return true;
        else
            return false;
    }
    else return true;  // can't enforce this for inter-chromosomals, so just return true
}


int getNonOverlap (const PAIR &a, const PAIR &b) {
    /* old, Contigger */
    // return ( abs(b.start1 - a.start1) + abs(b.start2 - a.start2) );
    return ( abs(b.core.start1 - a.core.start1) + abs(b.core.end2 - a.core.end2) );
}


bool doLengthsSupportOneAnother (const PAIR &a, const PAIR &b, int lengthDev)
{
    // cout << getFragSize(a) << "\t" << getFragSize(b) << "\t"
    //     << a.diffFromExpectedSize << "\t" << b.diffFromExpectedSize << "\t" << lengthDev << endl;
    return abs( getFragSize(a.core) - getFragSize(b.core) ) <= lengthDev;
}


//
bool doMultiLibraryLengthsSupportOneAnother (const PAIR &a, const PAIR &b, int aLengthDev, int bLengthDev)
{
    bool support = abs( (int) (a.diffFromExpectedSize - b.diffFromExpectedSize) ) <= max(aLengthDev, bLengthDev);

    // cout << getFragSize(a.core) << "\t" << getFragSize(b.core) << "\t"
    //      << a.diffFromExpectedSize << "\t" << b.diffFromExpectedSize << "\t"
    //      << aLengthDev << "\t" << bLengthDev << "\t" << support << endl;
    return support;
}


bool doSpansSupportOneAnother (const PAIR &a, const PAIR &b, int spanDev) {
    /********************************************************************************
      Two mapping spans support one another so long as there is at least 1bp of overlap
      __and__ the non-overlap between the spans is within the tolerance of the
      library(ies).  The latter constraint prevents situations where you have one
      base pair of overlap but the mappings are say, 100kb.  Exampli gratis:

      +...................................-
                                        +.........................................-

       [<-------------100kb----------------]
                                        [<-------------100kb----------------------]

       Yet, the two mappings come from a library where the median + 10(mad) is 2000bp.
       In such a case, the two mappings clearly do NOT support the same breakpoint.
    ***********************************************************************************/
    if ( (doSpansSupportCommonBreakpoint(a, b) == true) && (getNonOverlap(a, b) <= spanDev) ) {
        return true;
    }
    return false;
}


bool doMultiLibrarySpansSupportOneAnother (const PAIR &a, const PAIR &b, int aSpanDev, int bSpanDev) {
    /********************************************************************************
      Two mapping spans support one another so long as there is at least 1bp of overlap
      __and__ the non-overlap between the spans is within the tolerance of the
      library(ies).  The latter constraint prevents situations where you have one
      base pair of overlap but the mappings are say, 100kb.  Exampli gratis:

      +...................................-
                                        +.........................................-

       [<-------------100kb----------------]
                                        [<-------------100kb----------------------]

       Yet, the two mappings come from a library where the median + 10(mad) is 2000bp.
       In such a case, the two mappings clearly do NOT support the same breakpoint.
    ***********************************************************************************/
    if ( (doSpansSupportCommonBreakpoint(a, b) == true) &&
        (getNonOverlap(a, b) <= max(aSpanDev, bSpanDev)) )
    {
        return true;
    }
    return false;
}


int getTotalMM (const PAIR &pair) {
    return (pair.core.edit1 + pair.core.edit2);
}


int getTotalMMAmongAllMappings(const pairVector &mappings) {

    int totalMM = 0;

    // compute the total number of mismatches among all of the mappings
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter)
        totalMM += getTotalMM(*mapIter);

    return totalMM;
}


double getTotalWeightedSupportAmongAllMappings(const pairVector &mappings) {
    double totalWeightedSupport = 0.0;

    // compute the total weighted support among all of the mappings
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->mappingType == UNIQ_TYPE)
            totalWeightedSupport += UNIQ_WEIGHT;

        else if (mapIter->mappingType == ANCH_TYPE)
            totalWeightedSupport += ANCH_WEIGHT;

        else if (mapIter->mappingType == MULT_TYPE)
            totalWeightedSupport += MULT_WEIGHT;
    }
    return totalWeightedSupport;
}


// Return true if any of the mappings have include == true.
// That is, there is at least one mapping left in this putative cluster.
bool hasFinalSupport(const pairVector &mappings) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->include == true) return true;
    }
    return false;
}


int getClusterSizeAll(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->core.start1 < minStart1)    minStart1 = mapIter->core.start1;
        if (mapIter->core.end1 > maxEnd1)          maxEnd1 = mapIter->core.end1;
        if (mapIter->core.start2 < minStart2)    minStart2 = mapIter->core.start2;
        if (mapIter->core.end2 > maxEnd2)          maxEnd2 = mapIter->core.end2;
    }
    return (maxEnd2 - minStart1) + 1;
}


int getClusterSizeFinal(const pairVector &mappings, int &minStart1, int &minStart2, int &maxEnd1, int &maxEnd2) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        // only compute add to footprint if included in final call.
        if (mapIter->include == true) {
            if (mapIter->core.start1 < minStart1)    minStart1 = mapIter->core.start1;
            if (mapIter->core.end1 > maxEnd1)          maxEnd1 = mapIter->core.end1;
            if (mapIter->core.start2 < minStart2)    minStart2 = mapIter->core.start2;
            if (mapIter->core.end2 > maxEnd2)          maxEnd2 = mapIter->core.end2;
        }
    }
    return (maxEnd2 - minStart1) + 1;
}


bool isVariantUnlinked(const pairVector &mappings, int maxDistance) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if ( mapIter->core.chrom1 != mapIter->core.chrom2 ) return true;
        else if ( (mapIter->core.end2 - mapIter->core.start1) > maxDistance ) return true;
    }
    return false;
}


bool areAnyContigsUnlinked(const vector<IN_CLUSTER> &contigs) {
    vector<IN_CLUSTER>::const_iterator contigIter = contigs.begin();
    vector<IN_CLUSTER>::const_iterator contigEnd  = contigs.end();
    for (; contigIter != contigEnd; ++contigIter) {
        if ( contigIter->unlinked ) return true;
    }
    return false;
}


void getNumUniquePairs(const pairVector &mappings, int &numFinalUniques, int &numAllUniques) {

    // maps to track the unique pair ids that are clustered.
    std::map<string, short, std::less<string> > allUniqIds;
    std::map<string, short, std::less<string> > finalUniqIds;

    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->include == true) finalUniqIds[mapIter->core.readId]++;
        allUniqIds[mapIter->core.readId]++;
    }
    numFinalUniques = finalUniqIds.size();
    numAllUniques   = allUniqIds.size();
}


void getTotalEditDistance(const pairVector &mappings, int &totalMM1, int &totalMM2) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->include == true) {
            totalMM1 += mapIter->core.edit1;
            totalMM2 += mapIter->core.edit2;
        }
    }
}


void getTotalNumMappings(const pairVector &mappings, int &totalMappings1, int &totalMappings2) {
    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {
        if (mapIter->include == true) {
            if (mapIter->core.mate1 == 1) {
                totalMappings1 += mapIter->core.mappings1;
                totalMappings2 += mapIter->core.mappings2;
            }
            else {
                totalMappings1 += mapIter->core.mappings2;
                totalMappings2 += mapIter->core.mappings1;
            }
        }
    }
}


void computeSupport(const pairVector &mappings, int &finalSupport, double &finalWeightedSupport, double &allWeightedSupport,
                    int &numUniqueMappers, int &numAnchoredMappers, int &numMultipleMappers) {

    pairVector::const_iterator mapIter = mappings.begin();
    pairVector::const_iterator mapEnd  = mappings.end();
    for (; mapIter != mapEnd; ++mapIter) {

        if (mapIter->mappingType == UNIQ_TYPE) {
            if (mapIter->include) {
                finalSupport++;
                finalWeightedSupport += UNIQ_WEIGHT;
                numUniqueMappers++;
            }
            allWeightedSupport += UNIQ_WEIGHT;
        }
        else if (mapIter->mappingType == ANCH_TYPE) {
            if (mapIter->include) {
                finalSupport++;
                finalWeightedSupport += ANCH_WEIGHT;
                numAnchoredMappers++;
            }
            allWeightedSupport += ANCH_WEIGHT;
        }
        else if (mapIter->mappingType == MULT_TYPE) {
            if (mapIter->include) {
                finalSupport++;
                finalWeightedSupport += MULT_WEIGHT;
                numMultipleMappers++;
            }
            allWeightedSupport += MULT_WEIGHT;
        }
    }
}
