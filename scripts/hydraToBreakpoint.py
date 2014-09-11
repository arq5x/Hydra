#!/usr/bin/python

import sys
from optparse import OptionParser
import os

def hydraToBreakpoints(bedpe,slopOut,slopIn,localRange):
    #zero based/half-open
    slopIn += 1
    if bedpe == "stdin":
        hydraFile = sys.stdin
    else:
        hydraFile = open(str(bedpe))
    for line in hydraFile:
        split = line.split()
        chrom1 = split[0]
        start1 = int(split[1])
        end1 = int(split[2])
        chrom2 = split[3]
        start2 = int(split[4])
        end2 = int(split[5])
        ID = split[6]
        score = int(split[7])
        strand1 = split[8]
        strand2 = split[9]
        if chrom1 == chrom2 and start1 > start2:
            print "Error: malformed file " + str(start1) + " > " + str(start2) + " at id: " + ID
            sys.exit()
        if strand1 == strand2: #Classifying Inversions
            if chrom1  == chrom2 and end2 - start1 <= localRange:
                variant = "local_inversion"
            else:
                variant = "distant_inversion"
        elif strand1 == "+" and strand2 == "-": #Classifying Deletions
            if chrom1  == chrom2 and end2 - start1 <= localRange:
                variant = "local_deletion"
            else:
                variant = "distant_deletion"
        elif strand1 == "-" and strand2 == "+": #Classifying Tandem Duplications
            if chrom1  == chrom2 and end2 - start1 <= localRange:
                variant = "local_duplication"
            else:
                variant =  "distant_duplication"
        if strand1 == "+" and strand2 == "+":
            breakStart1 = str(end1-slopIn)
            breakEnd1 = str(end1+slopOut)
            breakStart2 = str(end2-slopIn)
            breakEnd2 = str(end2+slopOut)
        elif strand1 == "+" and strand2 == "-":
            breakStart1 = str(end1-slopIn)
            breakEnd1 = str(end1+slopOut)
            breakStart2 = str(start2-slopOut)
            breakEnd2 = str(start2+slopIn)
        elif strand1 == "-" and strand2 == "+":
            breakStart1 = str(start1-slopOut)
            breakEnd1 = str(start1+slopIn)
            breakStart2 = str(end2-slopIn)
            breakEnd2 = str(end2+slopOut)
        elif strand1 == "-" and strand2 == "-":
            breakStart1 = str(start1-slopOut)
            breakEnd1 = str(start1+slopIn)
            breakStart2 = str(start2-slopOut)
            breakEnd2 = str(start2+slopIn)
        print chrom1 + "\t" + breakStart1 + "\t" + breakEnd1 + "\t" + chrom2 + "\t" + breakStart2 + "\t" + breakEnd2 + "\t" + ID + "\t" + str(score) + "\t" + strand1 + "\t" + strand2 + "\t" + variant
    if bedpe != "stdin":
        hydraFile.close()
    
def main():
    usage = "%prog -b <bedpe> -c <contigFile.bedpe> [options]\nVersion: 0.1\nAuthor: Mitchell L. Leibowitz and Michael Linderg (Fixed Zero-based/half open)\n\n\
    Note: Add slopOut at your own risk.  If you add slopOut that goes across an entire deletion, you may get a false validation!"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="bedpe",
                        help="bedpe file in which entry 1 < entry 2; for stdin type \"stdin\"",
                        metavar="FILE")
    parser.add_option("-t", dest="slopOut", default=200, type="int", metavar="INT", help="slop out from the HYDRA footprint towards the breakpoint (half of average library size seems like a good number for predicting breakpoints) [default = 200]")
    parser.add_option("-n", dest="slopIn", default=10, type="int", metavar="INT", help="slop in towards the HYDRA footprint, away from the breakpoint [default = 10]")
    parser.add_option("-r", dest="localRange", default=1000000, type="int", metavar="INT", help="the range of bedpe coordinates considered local [default = 1000000]\nCalculated by subtracting field 6 from field 2.")
    (opts, args) = parser.parse_args()
    if opts.bedpe == None:
        parser.print_help()
        print
    else:
        hydraToBreakpoints(opts.bedpe, opts.slopOut, opts.slopIn, opts.localRange)

if __name__ == "__main__":
            sys.exit(main())
