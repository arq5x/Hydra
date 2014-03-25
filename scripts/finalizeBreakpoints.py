#!/usr/bin/env python
import sys
import os
import commands
from optparse import OptionParser


# constants
UNIQ_WEIGHT = 1.0
ANCH_WEIGHT = 0.5
MULT_WEIGHT = 0.1
UNIQ_TYPE   = 1
ANCH_TYPE   = 2
MULT_TYPE   = 3


class READ_TO_CLUSTER (object):
    """
    BEDPE class.  Can initialize to "NULL"
    by passign an empty list
    """ 
    def __init__(self, lineList = []):
        if len(lineList) > 0:
            self.chrom1                = lineList[0]
            self.start1                = int(lineList[1])
            self.end1                  = int(lineList[2])
            self.chrom2                = lineList[3]
            self.start2                = int(lineList[4])
            self.end2                  = int(lineList[5])
            self.readId                = lineList[6]
            self.whichMateIsBlock1     = int(lineList[7])
            self.strand1               = lineList[8]
            self.strand2               = lineList[9]
            self.edit1                 = int(lineList[10])
            self.edit2                 = int(lineList[11])
            self.mappings1             = int(lineList[12])
            self.mappings2             = int(lineList[13])
            self.mapq1                 = int(lineList[14])
            self.mapq2                 = int(lineList[15])
            self.mappingType           = int(lineList[16])
            self.diffFromExpectedSize  = int(lineList[17]) 
            self.support               = lineList[18]
            self.clusterId             = lineList[19]
            self.used                  = lineList[20]
            self.include               = lineList[21]
            self.fileNum               = lineList[22]
        else:
            self.valid = 0
            
    def report(self, out):
        out.write("\t".join([self.chrom1, str(self.start1), str(self.end1),
                             self.chrom2, str(self.start2), str(self.end2),
                             self.readId, str(self.whichMateIsBlock1), self.strand1, self.strand2,
                             str(self.edit1), str(self.edit2), str(self.mappings1), str(self.mappings2),
                             str(self.mapq1), str(self.mapq2), str(self.mappingType), str(self.diffFromExpectedSize),
                             self.support, self.clusterId, self.used, self.include, self.fileNum])
                )
        out.write("\n")
  
          
def getClusterSize(cluster):
    minStart1 = 999999999
    minStart2 = 999999999
    maxEnd1   = 0
    maxEnd2   = 0
    for r2c in cluster:
        if (r2c.start1 < minStart1):    minStart1 = r2c.start1
        if (r2c.end1 > maxEnd1):        maxEnd1   = r2c.end1
        if (r2c.start2 < minStart2):    minStart2 = r2c.start2
        if (r2c.end2 > maxEnd2):        maxEnd2   = r2c.end2
    return (maxEnd2 - minStart1) + 1


def getTotalMMAmongAllMappings(cluster):
    totalMM = 0
    for r2c in cluster:
        totalMM += r2c.edit1 + r2c.edit2
    return totalMM;


def getTotalWeightedSupportAmongAllMappings(cluster):
    totalWS = 0.0
    for r2c in cluster:
        if (r2c.mappingType == UNIQ_TYPE): 
            totalWS += UNIQ_WEIGHT
        elif (r2c.mappingType == ANCH_TYPE):
            totalWS += ANCH_WEIGHT
        elif (r2c.mappingType == MULT_TYPE):
            totalWS += MULT_WEIGHT
    return totalWS


def isVariantUnlinked(cluster, maxLinkedDistance):
    for r2c in cluster:
        if ( r2c.chrom1 != r2c.chrom2 ): 
            return 1
        elif ( (r2c.end2 - r2c.start1) > maxLinkedDistance ):
            return 1
    return 0


def getClusterSupport(cluster, maxLinkedDistance):
    contigSize         = getClusterSize(cluster)
    totalMM            = getTotalMMAmongAllMappings(cluster)
    weightedSupport    = getTotalWeightedSupportAmongAllMappings(cluster)
    unlinked           = isVariantUnlinked(cluster, maxLinkedDistance)
    return contigSize, totalMM, weightedSupport, unlinked
    
    
def computeSupportForEachCluster(clusterFile, maxLinkedDistance):
    """
    Rip through the master reads-to-clusters file
    and compute the support for each cluster. store
    this in a sqlite database.
    """

    prevClusterId = ""
    cluster = []
    clusterSupport = {}
    for line in open(clusterFile, 'r'):
        lineList = line.strip().split()
        r2c = READ_TO_CLUSTER(lineList)

        if (r2c.clusterId != prevClusterId and prevClusterId != ""):
            # grab the attributes 
            (size, totalMM, totalWS, unlinked) = getClusterSupport(cluster, maxLinkedDistance)
            clusterSupport[prevClusterId] = (prevClusterId, size, totalMM, totalWS, unlinked)
            cluster = []
            cluster.append(r2c)
        else:
            cluster.append(r2c)
        prevClusterId = r2c.clusterId
    
    (size, totalMM, totalWS, unlinked) = getClusterSupport(cluster, maxLinkedDistance)
    clusterSupport[prevClusterId] = (prevClusterId, size, totalMM, totalWS, unlinked)

    return clusterSupport
    

def sortClusterFileByReadId(clusterFile, sortMemory):
    outName = clusterFile + ".readsort"
    cmd = 'sort -S ' + sortMemory + ' -k7,7 ' + clusterFile + ' > ' + outName
    (status, output) = commands.getstatusoutput(cmd)
    return outName
    
def sortUpdatedFileByClusterId(clusterFile, sortMemory):        
    outName = clusterFile + ".clustersort"
    cmd = 'sort -S ' + sortMemory + ' -k20,20n ' + clusterFile + ' > ' + outName
    (status, output) = commands.getstatusoutput(cmd)
    return outName
    

def chooseBestClusterForRead(clusterStats):
    def bySupportThenMismatches(a,b):
        if a[3] > b[3]:
            return -1
        elif a[3] < b[3]:
            return 1
        elif a[2] < b[2]:
            return -1
        else:
            return 1
    if len(clusterStats) == 1:
        return clusterStats[0][0]
    else:
        distinct_support = list(set(clusterStats))
        distinct_support.sort(bySupportThenMismatches)
        return distinct_support[0][0]

def chooseBestClusterForReads(clusterFile, clusterSupport):
    """
    """
    def updateMappings(clusters, mappings, clusterSupport, out):
        support = []
        # collect the stats on all of the clusterIds that
        # this readpair belongs to
        for cluster in clusters:
            support.append(clusterSupport[cluster])
        # choose the best clusterId for this readpair.
        bestCluster = chooseBestClusterForRead(support)
        # update the membership of this read in all of the potential clusters.
        # That is, the best is included, the rest are excluded
        for mapping in mappings:
            if mapping.clusterId != bestCluster:
                mapping.include = 'N'
            else:
                mapping.include = 'Y'
            mapping.report(out)

    updatedClusterFile = clusterFile + ".updated"
    out = open(updatedClusterFile, 'w')
    
    prevReadId = ""
    clusters = []
    mappings = []

    for line in open(clusterFile, 'r'):
        lineList = line.strip().split()
        r2c = READ_TO_CLUSTER(lineList)

        # we've met a new read.  handle it if it belongs to
        # multiple clusters
        if (r2c.readId != prevReadId and prevReadId != ""):
            updateMappings(clusters, mappings, clusterSupport, out)
            # reset things for the next read
            clusters = []
            mappings = []
            clusters.append(r2c.clusterId)
            mappings.append(r2c)
        # same read.
        else:
            clusters.append(r2c.clusterId)
            mappings.append(r2c)
        prevReadId = r2c.readId
        
    # handle the last read in the file.
    updateMappings(clusters, mappings, clusterSupport, out)
        
    # We can also close the cursor if we are done with it
    out.close()
    return updatedClusterFile


def hasFinalSupport(mappings):
    for m in mappings:
        if m.include == "Y": 
            return True
    return False
    
def getClusterSizeFinal(mappings):
    minStart1            = minStart2 = 999999999
    maxEnd1              = maxEnd2   = 0
    for m in mappings:
        if (m.include == "Y"):
            if (m.start1 < minStart1):    minStart1 = m.start1
            if (m.end1 > maxEnd1):          maxEnd1 = m.end1
            if (m.start2 < minStart2):    minStart2 = m.start2
            if (m.end2 > maxEnd2):          maxEnd2 = m.end2
    return ((maxEnd2 - minStart1) + 1, minStart1, minStart2, maxEnd1, maxEnd2)

def getClusterSizeAll(mappings):
    minStart1            = minStart2 = 999999999
    maxEnd1              = maxEnd2   = 0
    for m in mappings:
        if (m.include == "Y"):
            if (m.start1 < minStart1):    minStart1 = m.start1
            if (m.end1 > maxEnd1):          maxEnd1 = m.end1
            if (m.start2 < minStart2):    minStart2 = m.start2
            if (m.end2 > maxEnd2):          maxEnd2 = m.end2
    return ((maxEnd2 - minStart1) + 1, minStart1, minStart2, maxEnd1, maxEnd2)

def getNumUniquePairs(mappings):
    # maps to track the unique pair ids that are clustered.
    allUniqIds   = {}
    finalUniqIds = {}

    for m in mappings:
        if (m.include == "Y"):
            if m.readId in finalUniqIds:
                finalUniqIds[m.readId] += 1
            else:
                finalUniqIds[m.readId] = 1
        
        if m.readId in allUniqIds:  
            allUniqIds[m.readId] += 1
        else:
            allUniqIds[m.readId] = 1
    numFinalUniques = len(finalUniqIds)
    numAllUniques   = len(allUniqIds)
    return (numFinalUniques, numAllUniques)


def getTotalEditDistance(mappings):
    totalMM1 = totalMM2 = 0
    for m in mappings:
        if (m.include == "Y"):
            totalMM1 += m.edit1
            totalMM2 += m.edit2
    return (totalMM1, totalMM2)
    
def getTotalMAPQ(mappings):
    totalMAPQ1 = totalMAPQ2 = 0
    for m in mappings:
        if (m.include == "Y"):
            totalMAPQ1 += m.mapq1
            totalMAPQ2 += m.mapq2
    return (totalMAPQ1, totalMAPQ2)

def getTotalNumMappings(mappings):
    totalMappings1 = totalMappings2 = 0
    for m in mappings:
        if (m.include == "Y"):
            totalMappings1 += m.mappings1
            totalMappings2 += m.mappings2
    return (totalMappings1, totalMappings2)
                
def computeSupport(mappings):
    allWeightedSupport   = 0.0
    finalSupport         = 0
    finalWeightedSupport = 0.0
    numUniqueMappers     = numAnchoredMappers = numMultipleMappers = 0
    
    for m in mappings:
        if (m.mappingType == UNIQ_TYPE):        
            if m.include == "Y": 
                finalSupport += 1
                finalWeightedSupport += UNIQ_WEIGHT
                numUniqueMappers += 1
            allWeightedSupport += UNIQ_WEIGHT
        elif (m.mappingType == ANCH_TYPE):
            if m.include == "Y":
                finalSupport += 1
                finalWeightedSupport += ANCH_WEIGHT
                numAnchoredMappers += 1
            allWeightedSupport += ANCH_WEIGHT
        elif (m.mappingType == MULT_TYPE):
            if m.include == "Y":
                finalSupport += 1
                finalWeightedSupport += MULT_WEIGHT
                numMultipleMappers += 1
            allWeightedSupport += MULT_WEIGHT
    return (finalSupport, finalWeightedSupport, allWeightedSupport, 
            numUniqueMappers, numAnchoredMappers, numMultipleMappers)
    
def reportCluster(mappings, contigId, final, all, detail):
    numMappings          = len(mappings)
    totalQual1           = totalQual2 = 0
    meanMappings1        = meanMappings2 = 0
    meanMAPQ1            = meanMAPQ2 = 0
    meanQual1            = meanQual2 = 0
    meanMM1              = meanMM2 = 0.0
    numAllUniques        = 0
    numFinalUniques      = 0
    contigSize           = 0
    
    # get the contig size and "footprints" based on the min and max mapping coordinates
    # if there is final support, only use the mappings that were _included_ in the cluster
    # otherwise, define the footprint based on all the mappings that were once in the cluster
    if (hasFinalSupport(mappings) == True):
        (contigSize, minStart1, minStart2, maxEnd1, maxEnd2) = getClusterSizeFinal(mappings)
    else:
        (contigSize, minStart1, minStart2, maxEnd1, maxEnd2)  = getClusterSizeAll(mappings)

    # (1)  Get _all_ of the unique ids in this cluster
    (numFinalUniques, numAllUniques) = getNumUniquePairs(mappings)
    
    # (2)  Get the total edit distance on each end of this cluster
    (totalMM1, totalMM2) = getTotalEditDistance(mappings)
    
    # (3)  Get the total number of mappings on each end of this cluster
    (totalMappings1, totalMappings2) = getTotalNumMappings(mappings)
    
    # (4)  Get the total MAPQ for each end of this cluster
    (totalMAPQ1, totalMAPQ2) = getTotalMAPQ(mappings)
    
    # (5)  Tally the final, finalWeighted and allWeightedSupport for this cluster
    (finalSupport, finalWeightedSupport, allWeightedSupport, 
     numUniqueMappers, numAnchoredMappers, numMultipleMappers) = computeSupport(mappings)

     
    # calculate quality statistics for each contig.
    if finalSupport > 0:
        meanMappings1 = totalMappings1 / finalSupport
        meanMappings2 = totalMappings2 / finalSupport
        meanQual1     = totalQual1 / finalSupport
        meanQual2     = totalQual2 / finalSupport
        meanMM1       = totalMM1 / finalSupport
        meanMM2       = totalMM2 / finalSupport
        meanMAPQ1     = totalMAPQ1 / finalSupport
        meanMAPQ2     = totalMAPQ2 / finalSupport
    
    # report this contigs as longs as there is at least one uniq
    # pair in the contig.
    if float(allWeightedSupport) > 0.0:

        contigId += 1;

        all.write("\t".join([mappings[0].chrom1, str(minStart1), str(maxEnd1), mappings[0].chrom2, str(contigId), 
                              str(finalSupport), mappings[0].strand1, mappings[0].strand2, str(meanMM1), str(meanMM2), 
                              str(meanMappings1), str(meanMappings2), str(meanMAPQ1), str(meanMAPQ2),
                              str(contigSize), str(numMappings), str(allWeightedSupport), 
                              str(finalSupport), str(finalWeightedSupport), str(numUniqueMappers), str(numAnchoredMappers), 
                              str(numMultipleMappers)]) + "\n")
        
        ### FIX: 2 should be replaced by a paramter indicating min support
        if (numFinalUniques >= 2):
            final.write("\t".join([mappings[0].chrom1, str(minStart1), str(maxEnd1), 
                                    mappings[0].chrom2, str(minStart2), str(maxEnd2),
                                    str(contigId), str(finalSupport), mappings[0].strand1, 
                                    mappings[0].strand2, str(meanMM1), str(meanMM2), 
                                    str(meanMappings1), str(meanMappings2), str(meanMAPQ1), str(meanMAPQ2),
                                    str(contigSize), str(numMappings), 
                                    str(allWeightedSupport), str(finalSupport), str(finalWeightedSupport), 
                                    str(numUniqueMappers), str(numAnchoredMappers), str(numMultipleMappers)]) + "\n")

            #write the read/mapping information for this contig to the contig detail file.
            for m in mappings:
                if (m.include == "Y"): 
                    detail.write("\t".join([m.chrom1, str(m.start1), str(m.end1), m.chrom2, 
                                             str(m.start2), str(m.end2), m.readId, str(m.whichMateIsBlock1), 
                                             m.strand1, m.strand2, str(m.edit1), str(m.edit2),
                                             str(m.mappings1), str(m.mappings2), str(m.mapq1), str(m.mapq2),
                                             str(m.mappingType), "Y", str(contigId)]) + "\n")
                else:
                    detail.write("\t".join([m.chrom1, str(m.start1), str(m.end1), m.chrom2, 
                                             str(m.start2), str(m.end2), m.readId, 
                                             str(m.whichMateIsBlock1), m.strand1, m.strand2, str(m.edit1), str(m.edit2),
                                             str(m.mappings1), str(m.mappings2), str(m.mapq1), str(m.mapq2),
                                             str(m.mappingType), "N", str(contigId)]) + "\n")
    return contigId


def createMasterAndDetailFiles(clusterSortedFile, outStub):
    
    finalFile  = outStub + ".final"
    allFile    = outStub + ".all"
    detailFile = outStub + ".detail"
    
    final    = open(finalFile, 'w')
    all      = open(allFile, 'w')
    detail   = open(detailFile, 'w')
    
    contigId = 0
    prevClusterId = ""
    mappings = []
        
    for line in open(clusterSortedFile, 'r'):
        lineList = line.strip().split()
        r2c = READ_TO_CLUSTER(lineList)

        # we've met a new read.  handle it if it belongs to
        # multiple clusters
        if (r2c.clusterId != prevClusterId and prevClusterId != ""):
            contigId = reportCluster(mappings, contigId, final, all, detail)
            
            # reset things for the next cluster
            mappings = []
            mappings.append(r2c)
        # same read.
        else:
            mappings.append(r2c)
        prevClusterId = r2c.clusterId
    # handle the last cluster in the file.
    contigId = reportCluster(mappings, contigId, final, all, detail)    

    final.close()
    all.close()
    detail.close()


 
def main():
    usage = """%prog -i <master>
    """
    parser = OptionParser(usage)
    
    parser.add_option("-i", dest="master", 
        help="Master assembled cluster file", 
        metavar="FILE")
    parser.add_option("-o", dest="outStub", 
        help="Stub for output files", 
        metavar="STRING")        
    parser.add_option("-m", action="store", type="string", dest="sortMemory", default="2G",
        help="The amount of memory to use for UNIX sort.  Default = 2G.  See \"man sort\" for valid values.", 
        metavar="STRING")
    parser.add_option("-d",action="store", type="int", dest="maxDist", default=1000000,
        help="Maximum intrachromosomal distance allowed before a variant is considered to be between unlinked DNA segments.", 
        metavar="INT")        
        
        
    # Grab the command line options
    (opts, args) = parser.parse_args()
    
    sortMemory = "2G"
    if opts.sortMemory is not None:
        sortMemory = opts.sortMemory
    
    # Make sure we have wha we need.  If so, on to the goodness.
    if opts.master is None:
        parser.print_help()
        print
    else:
        clusterSupport    = computeSupportForEachCluster(opts.master, opts.maxDist)
        readSortedFile    = sortClusterFileByReadId(opts.master, sortMemory)
        updatedFile       = chooseBestClusterForReads(readSortedFile, clusterSupport)
        clusterSortedFile = sortUpdatedFileByClusterId(updatedFile, sortMemory)
        createMasterAndDetailFiles(clusterSortedFile, opts.outStub)
        
        # clean up the temp files.
        cmd = 'rm ' + readSortedFile
        (status, output) = commands.getstatusoutput(cmd)
        cmd = 'rm ' + updatedFile
        (status, output) = commands.getstatusoutput(cmd)
        cmd = 'rm ' + clusterSortedFile
        (status, output) = commands.getstatusoutput(cmd)
            
if __name__ == "__main__":
    main()



