#!/usr/bin/env python
import sys
import commands
from optparse import OptionParser


class BEDPE_HYDRA (object):
    """
    BEDPE class.  Can initialize to "NULL"
    by passign an empty list
    """
    def __init__(self, bedList = []):
        if len(bedList) > 0:
            self.lineNum   = bedList[0]
            self.c1        = bedList[1]
            self.s1        = int(bedList[2])
            self.e1        = int(bedList[3])
            self.c2        = bedList[4]
            self.s2        = int(bedList[5])
            self.e2        = int(bedList[6])
            self.name      = bedList[7]
            self.score     = bedList[8]
            self.o1        = bedList[9]
            self.o2        = bedList[10]
            self.edit1     = int(bedList[11])
            self.edit2     = int(bedList[12])
            self.mappings1 = int(bedList[13])   # used to prioritize deduping.
            self.mappings2 = int(bedList[14])   # used to prioritize deduping.
            self.mapq1     = int(bedList[15])
            self.mapq2     = int(bedList[16])

            self.totalEdit     = self.edit1 + self.edit2
            self.totalMappings = self.mappings1 * self.mappings2
            self.valid = 1
        else:
            self.valid = 0

    def writeBed(self):
        print self.c1 + "\t" + str(self.s1) + "\t" + str(self.e1) + "\t" + \
            self.c2 + "\t" + str(self.s2) + "\t" + str(self.e2) + "\t" + \
            self.name + "\t" + self.score + "\t" + \
            self.o1 + "\t" + self.o2 + "\t" + \
            str(self.edit1) + "\t" + str(self.edit2) + "\t" + \
            str(self.mappings1) + "\t" + str(self.mappings2) + "\t" + \
            str(self.mapq1) + "\t" + str(self.mapq2)


def addLineNumAndSortFile(inbedpe, sortMemory, whichEnd):
    outFileName1 = inbedpe + ".sort1.lineNums.tmp"
    outFileName2 = inbedpe + ".sort2.lineNums.tmp"

    # Note that the columns being sorted are oner higher than the typical BEDPE columns
    # owing to the addition of the lineNumber as the first column via "cat -n"
    if whichEnd == 1:
        cmd = 'cat -n ' + inbedpe + ' | perl -lane \'if ($F[9] eq "-") {$temp = $F[3]; $F[3] = $F[2]; $F[2] = $temp;} if ($F[10] eq "-") {$temp = $F[6]; $F[6] = $F[5]; $F[5] = $temp;} print join("\t", @F);\' | sort -S  ' + sortMemory + ' -k2,2 -k5,5 -k3,3n -k6,6n -k10,10 -k11,11 > ' + outFileName1
        (status, output) = commands.getstatusoutput(cmd)
        return outFileName1
    else:
        cmd = 'sort -S ' + sortMemory + ' -k5,5 -k2,2 -k6,6n -k3,3n -k11,11 -k10,10 ' + outFileName1 + ' > ' + outFileName2
        (status, output) = commands.getstatusoutput(cmd)
        return outFileName2


def sortByName(inbedpe, sortMemory):
    outFileName = inbedpe + ".sortByName.lineNums.tmp"
    cmd = 'cat ' + inbedpe + ' | perl -lane \'if ($F[9] eq "-") {$temp = $F[3]; $F[3] = $F[2]; $F[2] = $temp;} if ($F[10] eq "-") {$temp = $F[6]; $F[6] = $F[5]; $F[5] = $temp;} print join("\t", @F);\' | sort -S  ' + sortMemory + ' -k8,8 > ' \
          + outFileName
    (status, output) = commands.getstatusoutput(cmd)
    return outFileName


def cleanUpFiles(fileList):
    for file in fileList:
        cmd = 'rm ' + file
        (status, output) = commands.getstatusoutput(cmd)


def byEditDistance(a, b):
    return cmp(a.totalEdit,b.totalEdit)


def isDuplicate(p, c, slop, needEither):
    """
    p = previous BEDPE entry.
    c = current BEDPE entry.

    1. are the mappings on the same chromosomes?
    2. are they the same strandedness?
    3. are they different read ids?
    4. are the positions within "slop" base pairs?

    If yes to all four, prev and curr are duplicates
    """
    if (needEither): # one or the other end must share a mapping
        if ((p.c1 == c.c1) and (p.c2 == c.c2) and (p.name != c.name) and
            (p.o1 == c.o1) and (p.o2 == c.o2) and
            ((abs(p.s1 - c.s1) <= slop) or (abs(p.s2 - c.s2) <= slop))):
            return 1
        else:
            return 0
    else: # both ends must share a mapping
        if ((p.c1 == c.c1) and (p.c2 == c.c2) and (p.name != c.name) and
            (p.o1 == c.o1) and (p.o2 == c.o2) and
            (abs(p.s1 - c.s1) <= slop) and (abs(p.s2 - c.s2) <= slop)):
            return 1
        else:
            return 0


def chooseBestMapping(dupList):
    """
    Sort the mappings in the duplicate block (asc) by total Edit distance.
    The one on top is the best, so return that.
    """
    dupList.sort(byEditDistance)
    return dupList[0]


def resolveDuplicateBlock(bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove, dupList, maxMappings):
    """
    Phase 1. Identify and remove all of the IDs having mappings greater than the allowable limit
    """
    for mapping in dupList:
        if mapping.totalMappings > maxMappings:
            idsToRemove[mapping.name] = True

    # remove the IDs with excessive mappings from further consideration
    dupList = [m for m in dupList if m.totalMappings <= maxMappings]


    """
    Phase 2. Now we need to determine whether or not an ID within this duplicate
             block has already been chosen as the BEST in a previously-examined block.
             If so, we will choose it again and make all of the other remaining IDs defer to it
             henceforth.
    """
    bestMapping = ''
    bestFound   = False
    for mapping in dupList:
        if mapping.name in bestDuplicates:
            bestMapping = mapping
            bestFound   = True
            break

    if bestFound == True:
        # Track the fact that these mappings should now defer to the best mapping.
        # This is used so that in future blocks, the other mappings from the inferior IDs
        # can defer to the best ID.
        for mapping in dupList:
            if mapping.name != bestMapping.name:
                inferiorDuplicates[mapping.name]  = bestMapping.name
                mappingsToRemove[mapping.lineNum] = True


    """
    Phase 3. By now, we know that a BEST mapping has NOT been previously found among the
             IDs in this block.  Now we will check if any of the IDs have been _replaced_
             by another better mapping in a previous block.
    """
    if bestFound == False:
        # sort the dup list be edit distance (asc) so that
        # the trumping mapping is found for IDs with the least edit distance first
        dupList.sort(byEditDistance)

        trumpMapping   = ''
        trumpedMapping = ''
        trumpFound   = False
        for mapping in dupList:
            if mapping.name in inferiorDuplicates:
                trumpMapping   = inferiorDuplicates[mapping.name]
                trumpedMapping = mapping.name
                trumpFound = True
                break
        if trumpFound == True:
            for mapping in dupList:
                if mapping.name != trumpedMapping:
                    inferiorDuplicates[mapping.name]  = trumpMapping
                    mappingsToRemove[mapping.lineNum] = True
        elif trumpFound == False and len(dupList) > 0:
            """
            Phase 4. If we have reached the phase, we know that we have encountered a block of IDs that
                     have never been seen before, nor have they been previously trumped by any other mapping.
                     In this scenario, we must choose the best mapping among the block and defer all of the
                     other mappings to it.
            """
            # By now, we know that the ids in this duplicate block are novel.
            # Thus, we choose the best among them and build a mapping between
            # it and the inferior mappings in this block
            bestMapping = chooseBestMapping(dupList)

            bestDuplicates[bestMapping.name] = True

            # Construct the mapping between the best mapping in this
            # block and all of the other ids that were not selected.
            for mapping in dupList:
                if mapping.name != bestMapping.name:
                    mappingsToRemove[mapping.lineNum] = True
                    if mapping.name not in inferiorDuplicates:
                        inferiorDuplicates[mapping.name] = bestMapping.name


def identifyDuplicates(bedpeIn, slop, either, maxMappings,
                       bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove):
    """
    Process the BEDPE file and flag duplicate ids.
    This function assumes that the BEDPE is sorted correctly.
    """
    currDupList         = []   # a list of current duplicate ids

    prev = BEDPE_HYDRA()    # create a null BEDPE for prev initially

    for line in open(bedpeIn, 'r'):
        lineList = line.strip().split()
        if (len(lineList) > 0):
            # create a BEDPE from the current line and if there
            # is a valid previous line, then check for duplicates
            curr = BEDPE_HYDRA(lineList)
            if (prev.valid == 1):
                if (isDuplicate(prev, curr, slop, either)):
                    currDupList.append(prev)
                    currDupList.append(curr)
                else:
                    # are there duplicates to resolve?
                    if (len(currDupList) > 0):
                        resolveDuplicateBlock(bestDuplicates, inferiorDuplicates, mappingsToRemove, \
                                              idsToRemove, currDupList, maxMappings)
                    # Reset the duplicate list for the next block
                    currDupList = []
            prev = curr

    # handle the last line
    if (len(currDupList) > 0):
        resolveDuplicateBlock(bestDuplicates, inferiorDuplicates, mappingsToRemove, \
                              idsToRemove, currDupList, maxMappings)

    return (bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove)


def reportCleanMappings(bedpeIn, inferiorDuplicates,
                        mappingsToRemove, idsToRemove):
    for line in open(bedpeIn, 'r'):
        lineList = line.strip().split()
        curr = BEDPE_HYDRA(lineList)
        if (curr.name    not in idsToRemove and  \
            curr.lineNum not in mappingsToRemove):
            if (curr.name not in inferiorDuplicates):
                curr.writeBed()
            else:
                curr.name = inferiorDuplicates[curr.name]
                curr.writeBed()

def main():
    usage = """%prog -i <BEDPE> -s <INT, default=1>
    """
    parser = OptionParser(usage)

    parser.add_option("-i", dest="inbedpe",
        help="BEDPE input file",
        metavar="FILE")

    parser.add_option("-s",action="store", type="int", dest="slop", default=1,
        help="The amount of slop, in bp, allowed for duplicates to be marked",
        metavar="INT")

    parser.add_option("-m", action="store", type="int", dest="maxMappings", default=1,
        help="The maximum number of mapping combinations allowed before a duplicate ID has _all_ of its mappings removed",
        metavar="INT")

    parser.add_option("--mem", action="store", type="string", dest="sortMemory", default="2G",
        help="The amount of memory to use for UNIX sort.  Default = 2G.  See \"man sort\" for valid values.",
        metavar="STRING")

    parser.add_option("-e", dest="either", action="store_true",
        help="Dedup if __either__ end shares a start position (within slop) in common.")

    # Grab the command line options
    (opts, args) = parser.parse_args()

    sortMemory = "12G"
    if opts.sortMemory is not None:
        sortMemory = opts.sortMemory

    # Make sure we have wha we need.  If so, on to the goodness.
    if opts.inbedpe is None:
        parser.print_help()
        print
    else:

        bestDuplicates      = {}  # mapping of the best id chosen in a block
        inferiorDuplicates  = {}  # mapping of inferior ids to the best mapping chosen
        mappingsToRemove    = {}  # lineNums that were in duplicate blocks and were inferior.
        idsToRemove         = {}  # ids that were in duplicate blocks and had excessive mappings.
                                  # remove of their mappings.

        tmpFile1 = addLineNumAndSortFile(opts.inbedpe, sortMemory, 1)

        # find duplicates using a sort focused on end1
        (bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove) = \
             identifyDuplicates(tmpFile1, opts.slop, opts.either, opts.maxMappings,
                                bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove)

        if (opts.either == True):
            tmpFile2 = addLineNumAndSortFile(opts.inbedpe, sortMemory, 2)
            (bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove) = \
                 identifyDuplicates(tmpFile2, opts.slop, opts.either, opts.maxMappings,
                                    bestDuplicates, inferiorDuplicates, mappingsToRemove, idsToRemove)

            tmpFileSortedByName = sortByName(tmpFile2, sortMemory)
            reportCleanMappings(tmpFileSortedByName, inferiorDuplicates, mappingsToRemove, idsToRemove)

            # delete the temp files
            cleanUpFiles([tmpFile1, tmpFile2, tmpFileSortedByName])
        else:
            tmpFileSortedByName = sortByName(tmpFile1, sortMemory)
            reportCleanMappings(tmpFileSortedByName, inferiorDuplicates, mappingsToRemove, idsToRemove)
            # delete the temp files
            cleanUpFiles([tmpFile1, tmpFileSortedByName])


if __name__ == "__main__":
    main()



