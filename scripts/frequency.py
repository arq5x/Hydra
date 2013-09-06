#!/usr/local/bin/python2.6                                                                                                                                                                                               

import sys
import collections
from optparse import OptionParser

def mapSamples(samplesFile):
    """
        Create a map of all the sample names and their 
        0-based position in the samples file.  This position
        will be used to drive the location of counts for this sample
        in a list used to store counts for all N samples.
    """
    sampleMap = collections.defaultdict(int)
    sampleList = []
    for i, line in enumerate(open(samplesFile, 'r')):
	sample = line.strip().split()[0]
	sampleList.append(sample)
        sampleMap[sample] = i
    return sampleMap, sampleList


def printHeader(sampleMap, sampleList):
    headerRoot = ["#chrom1_1", "minStart1_2", "maxEnd1_3", 
                  "chrom2_4", "minStart2_5", "maxEnd2_6", 
                  "contigId_7", "distinctPairs_8", "strand1_9", 
                  "strand2_10", "meanMM1_11", "meanMM2_12", 
                  "meanMappings1_13", "meanMappings2_14", 
                  "meanMAPQ1_15", "meanMAPQ2_16", "contigSize_17" ,
                  "numPairs_18", "allWeightedSupport_19", "finalSupport_20", 
                  "finalWeightedSupport_21", "numUniq_22", "numArch_23", "numMult_24"]
    print "\t".join(headerRoot),
    print "\t",
    print "Y\t".join(sampleList) + "Y"
    
    
def countSamplesPerBreakpoint(detailFile, finalFile, sampleMap, sampleList):
    
    def count(idx, yn, countsY):
        if (yn == "Y"):
            counts_Y[idx]  += 1
    

    def getSample(readId, samples):
    	for i in samples:
	    s = readId.find(i)
	    if s == 0:
	       e = len(i)
	       return i[s:e]

    numSamples = len(sampleMap)
    counts_Y = [0]*numSamples        
    prevBreak = 'init'
    final = open(finalFile, 'r')
    for line in open(detailFile, 'r'):
        det = line.strip().split()
        sample = getSample(det[6], sampleList)
        idx = sampleMap[sample]
        yn = det[17] # ARQ s/15/17/, to reflect new cols. in detail
        currBreak = det[18] # ARQ s/16/18/, to reflect new cols. in detail
        if (currBreak == prevBreak) or (prevBreak == 'init'):
            count(idx, yn, counts_Y)
        else:
            finalLine = final.readline()
            reportFrequencies(finalLine, counts_Y, prevBreak)
            counts_Y = [0]*numSamples
            count(idx, yn, counts_Y)
        prevBreak = currBreak
    finalLine = final.readline()
    reportFrequencies(finalLine, counts_Y, prevBreak)


def reportFrequencies(finalLine, counts_Y, prevBreak):
    finalLineList = finalLine.strip().split()
    print finalLine.strip() + "\t" + "\t".join("%s" % c for c in counts_Y)
    
            
def main():
    parser = OptionParser()
    parser.add_option("-d", dest="detail", help="detail filename", metavar="FILE")
    parser.add_option("-f", dest="final", help="final filename", metavar="FILE")
    parser.add_option("-s", dest="samples", help="samples filename", metavar="FILE")
    parser.add_option("-x", action="store_false", dest="col", default=True, help="do not print column headers")
    (options, args) = parser.parse_args()
    # build a map of the samples
    (sampleMap, sampleList) = mapSamples(options.samples)
    printHeader(sampleMap, sampleList)
    countSamplesPerBreakpoint(options.detail, options.final, sampleMap, sampleList)

        
        
if __name__ == "__main__":
        sys.exit(main())
