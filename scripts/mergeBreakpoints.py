#!/usr/bin/env python
import sys
import os
import subprocess
import commands
import copy
import shutil
from optparse import OptionParser


class HYDRA_FINAL (object):
    def __init__(self, lineList = []):
        if (len(lineList) > 0):
            self.chrom1                  = lineList[0]
            self.start1                  = int(lineList[1])
            self.end1                    = int(lineList[2])
            self.chrom2                  = lineList[3]
            self.start2                  = int(lineList[4])
            self.end2 	                 = int(lineList[5])
            self.breakpointId 	         = lineList[6]
            self.numDistinctPairs 	     = int(lineList[7])
            self.strand1                 = lineList[8]
            self.strand2                 = lineList[9]
            self.meanEditDist1           = int(lineList[10])
            self.meanEditDist2           = int(lineList[11])
            self.meanMappings1           = int(lineList[12])
            self.meanMappings2           = int(lineList[13])
            self.breakpointSize          = int(lineList[14])
            self.numMappings             = int(lineList[15])
            self.allWeightedSupport 	 = float(lineList[16])
            self.finalSupport            = int(lineList[17])
            self.finalWeightedSupport    = float(lineList[18])
            self.numUniquePairs          = int(lineList[19])
            self.numAnchoredPairs 	     = int(lineList[20])
            self.numMultiplyMappedPairs  = int(lineList[21])
            self.valid = True
        else:
            self.valid = False

    def __eq__(self, other):
        return self.breakpointId == other.breakpointId

    def dump(self, out):
        if self.valid:
            out.write("\t".join([self.chrom1, str(self.start1), str(self.end1),
                             self.chrom2, str(self.start2), str(self.end2),
                             str(self.breakpointId), str(self.numDistinctPairs), self.strand1, self.strand2,
                             str(self.meanEditDist1), str(self.meanEditDist2), str(self.meanMappings1), str(self.meanMappings2),
                             str(self.breakpointSize), str(self.numMappings), str(self.allWeightedSupport), str(self.finalSupport),
                             str(self.finalWeightedSupport), str(self.numUniquePairs), str(self.numAnchoredPairs), str(self.numMultiplyMappedPairs)]))
            out.write("\n")
            
class HYDRA_FINAL_OVERLAP (object):
    def __init__(self, lineList = []):
        if (len(lineList) > 0):
            self.a = HYDRA_FINAL(lineList[0:22])
            self.b = HYDRA_FINAL(lineList[22:44])
            
            self.a_chrom1                  = lineList[0]
            self.a_start1                  = int(lineList[1])
            self.a_end1                    = int(lineList[2])
            self.a_chrom2                  = lineList[3]
            self.a_start2                  = int(lineList[4])
            self.a_end2                  = int(lineList[5])
            self.a_breakpointId          = lineList[6]
            self.a_numDistinctPairs      = int(lineList[7])
            self.a_strand1                 = lineList[8]
            self.a_strand2                 = lineList[9]
            self.a_meanEditDist1           = int(lineList[10])
            self.a_meanEditDist2           = int(lineList[11])
            self.a_meanMappings1           = int(lineList[12])
            self.a_meanMappings2           = int(lineList[13])
            self.a_breakpointSize          = int(lineList[14])
            self.a_numMappings             = int(lineList[15])
            self.a_allWeightedSupport        = float(lineList[16])
            self.a_finalSupport            = int(lineList[17])
            self.a_finalWeightedSupport    = float(lineList[18])
            self.a_numUniquePairs          = int(lineList[19])
            self.a_numAnchoredPairs      = int(lineList[20])
            self.a_numMultiplyMappedPairs  = int(lineList[21])
            
            self.b_chrom1                  = lineList[22]
            self.b_start1                  = int(lineList[23])
            self.b_end1                    = int(lineList[24])
            self.b_chrom2                  = lineList[25]
            self.b_start2                  = int(lineList[26])
            self.b_end2                  = int(lineList[27])
            self.b_breakpointId          = lineList[28]
            self.b_numDistinctPairs      = int(lineList[29])
            self.b_strand1                 = lineList[30]
            self.b_strand2                 = lineList[31]
            self.b_meanEditDist1           = int(lineList[32])
            self.b_meanEditDist2           = int(lineList[33])
            self.b_meanMappings1           = int(lineList[34])
            self.b_meanMappings2           = int(lineList[35])
            self.b_breakpointSize          = int(lineList[36])
            self.b_numMappings             = int(lineList[37])
            self.b_allWeightedSupport        = float(lineList[38])
            self.b_finalSupport            = int(lineList[39])
            self.b_finalWeightedSupport    = float(lineList[40])
            self.b_numUniquePairs          = int(lineList[41])
            self.b_numAnchoredPairs      = int(lineList[42])
            self.b_numMultiplyMappedPairs  = int(lineList[43])

            self.valid = True
        else:
            self.a = HYDRA_FINAL([])
            self.b = HYDRA_FINAL([])
            self.a_breakpointId = -1
            self.valid = False

    def dumpA(self, out):
        if self.valid:
            out.write("\t".join([self.a.chrom1, str(self.a.start1), str(self.a.end1),
                         self.a.chrom2, str(self.a.start2), str(self.a.end2),
                         str(self.a.breakpointId), str(self.a.numDistinctPairs), self.a.strand1, self.a.strand2,
                         str(self.a.meanEditDist1), str(self.a.meanEditDist2), str(self.a.meanMappings1), str(self.a.meanMappings2),
                         str(self.a.breakpointSize), str(self.a.numMappings), str(self.a.allWeightedSupport), str(self.a.finalSupport),
                         str(self.a.finalWeightedSupport), str(self.a.numUniquePairs), str(self.a.numAnchoredPairs), str(self.a.numMultiplyMappedPairs)]))
            out.write("\n")

    def dumpB(self, out):
        if self.valid:
            out.write("\t".join([self.b.chrom1, str(self.b.start1), str(self.b.end1),
                         self.b.chrom2, str(self.b.start2), str(self.b.end2),
                         str(self.b.breakpointId), str(self.b.numDistinctPairs), self.b.strand1, self.b.strand2,
                         str(self.b.meanEditDist1), str(self.b.meanEditDist2), str(self.b.meanMappings1), str(self.b.meanMappings2),
                         str(self.b.breakpointSize), str(self.b.numMappings), str(self.b.allWeightedSupport), str(self.b.finalSupport),
                         str(self.b.finalWeightedSupport), str(self.b.numUniquePairs), str(self.b.numAnchoredPairs), str(self.b.numMultiplyMappedPairs)]))
            out.write("\n")


    def debug(self, out):
        if self.valid:
            out.write("\t".join([self.a.chrom1,
                                 str(self.a.start1),
                                 str(self.a.end1),
                                 self.a.chrom2,
                                 str(self.a.start2),
                                 str(self.a.end2),
                                 str(self.a.breakpointId),
                                 str(self.a.numDistinctPairs),
                                 self.a.strand1,
                                 self.a.strand2,
                         self.b.chrom1, str(self.b.start1), str(self.b.end1),
                         self.b.chrom2, str(self.b.start2), str(self.b.end2),
                         str(self.b.breakpointId), str(self.b.numDistinctPairs), self.b.strand1, self.b.strand2]))





def overlap(s1, s2, e1, e2):
    return min(e1,e2) - max(s1,s2)

def innerSpanOverlap(p,c):
    if p.valid == False:
        return -1
    else:
        return overlap(p.end1, c.end1, p.start2, c.start2)

def doFootPrintsSupport(p,c,slop):
    # both footprints overlap
    left  = overlap(p.start1, c.start1, p.end1, c.end1)
    right = overlap(p.start2, c.start2, p.end2, c.end2)
    if (left > (-1 * slop)) and (right > (-1 * slop)):
        return True
    else:
        return False

def supportSameBreakpoint(a, b, slop):
    if (a.valid == False) or (b.valid == False):
        return False
    if ((a.chrom1 == a.chrom2) and (b.chrom1 == b.chrom2)):
        return (innerSpanOverlap(a,b) > 0) and (doFootPrintsSupport(a,b,slop) == True)
    else:
        return doFootPrintsSupport(a,b,slop)


def merge_s(breaks, output, slop):
    # Step 1. Ensure that all breaks in the set support the same breakpoint.
    toMerge = []
    for b1 in breaks:
        allInSupport = True
        for b2 in breaks:
            if (supportSameBreakpoint(b1, b2, slop) == False):
                allInSupport = False
        if allInSupport == True:
            toMerge.append(b1)
        else:
            b1.dump(output)

    if len(toMerge) <= 1:
        for b in toMerge:
            b.dump(output)
        return False

    minStart1 = minStart2      = 999999999
    maxEnd1   = maxEnd2        = 0
    ids                        = []
    maxId                      = 0
    allNumDistinctPairs 	   = 0
    allMeanEditDist1           = 0
    allMeanEditDist2           = 0
    allMeanMappings1           = 0
    allMeanMappings2           = 0
    allNumMappings             = 0
    allAllWeightedSupport 	   = 0
    allFinalSupport            = 0
    allFinalWeightedSupport    = 0
    allNumUniquePairs          = 0
    allNumAnchoredPairs 	   = 0
    allNumMultiplyMappedPairs  = 0

    numBreaks = len(toMerge)
    for b in toMerge:
        # establish the coordinates
        if (b.start1 < minStart1):    minStart1 = b.start1
        if (b.end1 > maxEnd1):          maxEnd1 = b.end1
        if (b.start2 < minStart2):    minStart2 = b.start2
        if (b.end2 > maxEnd2):          maxEnd2 = b.end2
        ids.append(str(b.breakpointId))
        allNumDistinctPairs 	   += b.numDistinctPairs
        allMeanEditDist1           += b.meanEditDist1
        allMeanEditDist2           += b.meanEditDist2
        allMeanMappings1           += b.meanMappings1
        allMeanMappings2           += b.meanMappings2
        allNumMappings             += b.numMappings
        allAllWeightedSupport 	   += b.allWeightedSupport
        allFinalSupport            += b.finalSupport
        allFinalWeightedSupport    += b.finalWeightedSupport
        allNumUniquePairs          += b.numUniquePairs
        allNumAnchoredPairs 	   += b.numAnchoredPairs
        allNumMultiplyMappedPairs  += b.numMultiplyMappedPairs

    finalMeanEditDist1 = allMeanEditDist1 / numBreaks
    finalMeanEditDist2 = allMeanEditDist2 / numBreaks
    finalMeanMappings1 = allMeanMappings1 / numBreaks
    finalMeanMappings2 = allMeanMappings2 / numBreaks
    mergedSize = (maxEnd2 - minStart1) + 1

    idString = "_".join(ids)
    mergedBreak = HYDRA_FINAL([breaks[0].chrom1, str(minStart1), str(maxEnd1),
                               breaks[0].chrom2, str(minStart2), str(maxEnd2),
                               idString, str(allNumDistinctPairs), breaks[0].strand1, breaks[0].strand2,
                               str(finalMeanEditDist1), str(finalMeanEditDist2), str(finalMeanMappings1), str(finalMeanMappings2),
                               str(mergedSize), str(allNumMappings), str(allAllWeightedSupport), str(allFinalSupport),
                               str(allFinalWeightedSupport), str(allNumUniquePairs), str(allNumAnchoredPairs), str(allNumMultiplyMappedPairs)])
    mergedBreak.dump(output)
    return True



def merge(toMerge, output):
    # First, capture the info from the "A" breakpoint.
    minStart1 = minStart2      = 999999999
    maxEnd1   = maxEnd2        = 0
    ids                        = []
    maxId                      = 0
    allNumDistinctPairs 	   = 0
    allMeanEditDist1           = 0
    allMeanEditDist2           = 0
    allMeanMappings1           = 0
    allMeanMappings2           = 0
    allNumMappings             = 0
    allAllWeightedSupport 	   = 0
    allFinalSupport            = 0
    allFinalWeightedSupport    = 0
    allNumUniquePairs          = 0
    allNumAnchoredPairs 	   = 0
    allNumMultiplyMappedPairs  = 0

    numBreaks = len(toMerge)

    for b in toMerge:

        # establish the coordinates
        if (b.start1 < minStart1):    minStart1 = b.start1
        if (b.end1 > maxEnd1):          maxEnd1 = b.end1
        if (b.start2 < minStart2):    minStart2 = b.start2
        if (b.end2 > maxEnd2):          maxEnd2 = b.end2
        ids.append(str(b.breakpointId))

        allNumDistinctPairs 	   += b.numDistinctPairs
        allMeanEditDist1           += b.meanEditDist1
        allMeanEditDist2           += b.meanEditDist2
        allMeanMappings1           += b.meanMappings1
        allMeanMappings2           += b.meanMappings2
        allNumMappings             += b.numMappings
        allAllWeightedSupport 	   += b.allWeightedSupport
        allFinalSupport            += b.finalSupport
        allFinalWeightedSupport    += b.finalWeightedSupport
        allNumUniquePairs          += b.numUniquePairs
        allNumAnchoredPairs 	   += b.numAnchoredPairs
        allNumMultiplyMappedPairs  += b.numMultiplyMappedPairs

    finalMeanEditDist1 = allMeanEditDist1 / numBreaks
    finalMeanEditDist2 = allMeanEditDist2 / numBreaks
    finalMeanMappings1 = allMeanMappings1 / numBreaks
    finalMeanMappings2 = allMeanMappings2 / numBreaks
    mergedSize = (maxEnd2 - minStart1) + 1

    idString = "_".join(ids)
    mergedBreak = HYDRA_FINAL([toMerge[0].chrom1, str(minStart1), str(maxEnd1),
                               toMerge[0].chrom2, str(minStart2), str(maxEnd2),
                               idString, str(allNumDistinctPairs), toMerge[0].strand1, toMerge[0].strand2,
                               str(finalMeanEditDist1), str(finalMeanEditDist2), str(finalMeanMappings1), str(finalMeanMappings2),
                               str(mergedSize), str(allNumMappings), str(allAllWeightedSupport), str(allFinalSupport),
                               str(allFinalWeightedSupport), str(allNumUniquePairs), str(allNumAnchoredPairs), str(allNumMultiplyMappedPairs)])
    mergedBreak.dump(output)
    
    
def merge_p2p(breaks, output, slop, breaksSeen, seenButNotWritten):

    mergingDone = False
    for b1 in breaks:
        # skip if we have already dealt with this breakpoint
        if b1.breakpointId not in breaksSeen:
            for b2 in breaks:
                if b1.breakpointId != b2.breakpointId and b2.breakpointId not in breaksSeen:
                    if (supportSameBreakpoint(b1, b2, slop) == True):
                        breaksSeen[b1.breakpointId] = 1
                        breaksSeen[b2.breakpointId] = 1
                        merge([b1, b2], output)
                        mergingDone = True
                        break
        if mergingDone:
            break

    for b1 in breaks:
        if b1.breakpointId not in breaksSeen:
            seenButNotWritten.append(b1)

    return mergingDone


def mergeBreakpointsWithSorting(finalFile, slop, sortMethod):


    prev = HYDRA_FINAL([])
    set  = []

    cmd = 'LC_ALL=C sort -k9,9 -k10,10 -k2,2n ' + finalFile # default to sort by start (assumes already sorted by chr/chr)
    if (sortMethod == "sortByEnd"):
        cmd = 'LC_ALL=C sort -k9,9 -k10,10 -k5,5n ' + finalFile # default to sort by start (assumes already sorted by chr/chr)

    sortedFile = os.popen(cmd)
    somethingMerged = False

    # make a new final file.
    outputFile = "." + finalFile + "." + sortMethod
    output = open(outputFile, 'w')
    prevOverlap = False
    for line in sortedFile:
        lineList = line.strip().split("\t")
        curr = HYDRA_FINAL(lineList)

        if (supportSameBreakpoint(prev, curr, slop)):
            set.append(prev)
            prevOverlap = True
        else:
            if prevOverlap == True:
                set.append(prev)
                merged = merge_s(set, output, slop)
                if (merged == True):
                    somethingMerged = True
                prevOverlap = False
                set = []
            else:
                prev.dump(output)
                set = []
        prev = curr

    # handle last set
    if prevOverlap == True:
        set.append(curr)
    if len(set) > 1:
        merge_s(set, output, slop)
    else:
        prev.dump(output)
            
    output.close()
    return (outputFile, somethingMerged)



def mergeBreakpointsWithPairToPair(finalFile, slop, round):

    breaksSeen = {}
    seenButNotWritten = []
    prev = HYDRA_FINAL_OVERLAP([])
    set  = []
    cmd = 'pairToPair -rdn -a ' + finalFile + ' -b ' + finalFile + ' -slop ' + str(slop)
    p2p = os.popen(cmd)
    # make a new final file.
    outputFile = "." + finalFile + ".p2p"
    output = open(outputFile, 'w')

    somethingMerged = False
    for line in p2p:
        lineList = line.strip().split("\t")
        curr = HYDRA_FINAL_OVERLAP(lineList)
        if (curr.a_breakpointId != prev.a_breakpointId):
            if len(set) > 0:
                breakList = []
                for br in set:
                    breakList.append(br.a)
                    breakList.append(br.b)

                merged = merge_p2p(breakList, output, slop, breaksSeen, seenButNotWritten)
                if merged == True:
                    somethingMerged = True
                set = []
        set.append(curr)
        prev = curr

    # handle last set in the pairs that have overlap
    breakList = []
    for br in set:
        breakList.append(br.a)
        breakList.append(br.b)
    merged = merge_p2p(breakList, output, slop, breaksSeen, seenButNotWritten)
    if merged == True:
        somethingMerged = True

    # handle the breakpoints that were caught bu p2p, but did not support the same breakpoint.
    # write out the ones that have not yet been written
    for b in seenButNotWritten:
        if b.breakpointId not in breaksSeen:
            b.dump(output)
            breaksSeen[b.breakpointId] = 1

    # report the pairs that had not overlap
    cmd = 'pairToPair -a ' + finalFile + ' -b ' + finalFile + ' -slop ' + str(slop) + ' -rdn -type notboth'
    p2p = os.popen(cmd)
    for line in p2p:
        lineList = line.strip().split("\t")
        curr = HYDRA_FINAL(lineList)
        curr.dump(output)

    # close up shop
    output.close()
    return (outputFile, somethingMerged)





def main():
    usage = """%prog -i <master>
    """
    parser = OptionParser(usage)

    parser.add_option("-f", dest="final",
        help="Hydra .final file",
        metavar="FILE")
    parser.add_option("-d", dest="detail",
        help="Hydra .detail file",
        metavar="FILE")
    parser.add_option("-o", dest="outStub",
        help="Merged final file",
        metavar="STRING")
    parser.add_option("-s",action="store", type="int", dest="maxDist", default=100,
        help="Maximum distance allowed for two flanking breakpoints to be merged.",
        metavar="INT")


    # Grab the command line options
    (opts, args) = parser.parse_args()

    # Make sure we have wha we need.  If so, on to the goodness.
    if opts.final is None:
        parser.print_help()
        print
    else:

        round = 0
        mergingDone = True
        finalInput = opts.final
        mergeByStartAndEnd = None
        tempFiles = []

        # Phase 1. Merge by sorting
        while (mergingDone == True and round <= 20):
            (mergeByStart, mergedByStart)       = mergeBreakpointsWithSorting(finalInput,  opts.maxDist,  "sortByStart")
            (mergeByStartAndEnd, mergedByEnd)   = mergeBreakpointsWithSorting(mergeByStart, opts.maxDist, "sortByEnd")
            # Was merging done in either case?  If so, we want to resort and try again.
            mergingDone = (mergedByStart or mergedByEnd)
            if mergingDone:
                round += 1
                finalInput = "." + opts.final + "." +  str(round)
                shutil.copy(mergeByStartAndEnd, finalInput)
            tempFiles.append(finalInput)
            tempFiles.append(mergeByStart)
            tempFiles.append(mergeByStartAndEnd)

        # Phase 2. Merge by pairing
        mergingDone = True
        mergeByP2P = None
        finalInput = mergeByStartAndEnd
        while (mergingDone == True):
            (mergeByP2P, mergingDone) = mergeBreakpointsWithPairToPair(finalInput,  opts.maxDist, round)
            if mergingDone:
                round += 1
                finalInput = "." + opts.final + "." +  str(round)
                shutil.copy(mergeByP2P, finalInput)
            tempFiles.append(finalInput)
            tempFiles.append(mergeByP2P)

        # Phase 3. Remap the final and detail files
        mergedFinal  = opts.outStub + ".final"
        mergeMapping = opts.outStub + ".mergemap"
        final = open(mergedFinal, 'w')
        mmap  = open(mergeMapping, 'w')

        # map of old ids to new, merged ids
        oldToNew = {}
        # uniqify the breakpoints:
        for line in open(finalInput, 'r'):
            lineList = line.strip().split()
            ids = lineList[6].split("_")
            uniqIds = set(ids)
            # just mark a merged breakpoint
            uniqMaster  = list(uniqIds)[0]
            if len(uniqIds) > 1:
                uniqMaster += ".m"
            # replace the long list of merged IDs with a single, representative ID
            lineList[6] = uniqMaster
            # save the mappings of individual ids to the repr. ID
            mmap.write(uniqMaster + "\t" + ",".join(list(uniqIds)) + "\n")

            for id in uniqIds:
                oldToNew[id] = uniqMaster
            final.write("\t".join(lineList))
            final.write("\n")
        final.close()
        mmap.close()

        # now the detail file.
        mergedDetail = opts.outStub + ".detail"
        detail = open(mergedDetail, 'w')
        for line in open(opts.detail, 'r'):
            lineList = line.strip().split("\t")

            # only replace the id if it ever made it to the ".final" file.
            # i.e., we don't care about "N" mappings, just "Y"
            # NOTE: The "if oldId in oldToNew" check is a hack to correct for the
            # fact that I mistakenly wrote breakpoints with insufficient final support
            # to the detail file.  This check will correct this problem.
            oldId = lineList[16]
            newId = oldToNew[oldId]
            lineList[16] = newId
            detail.write("\t".join(lineList))
            detail.write("\n")
            # if oldId in oldToNew:
            #     
            #     if lineList[15] == "Y":
            #         newId = oldToNew[oldId]
            #         lineList[16] = newId
            #         detail.write("\t".join(lineList))
            #     else:
            #         detail.write("\t".join(lineList))
            #     detail.write("\n")
            # else:
            #     print oldId, "not in set"
        detail.close()
        
        for tmp in tempFiles:
            # cleanup temp files and exit
            if os.path.exists(tmp):
                os.remove(tmp)
        
        print "Finished."

if __name__ == "__main__":
    main()



