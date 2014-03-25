#!/usr/bin/env python
import sys
import getopt
import string
import os
from optparse import OptionParser
import shutil

def processBEDPE(bedpe, dist, out, genome):
    
    bed12 = bedpe + ".bed12"
    cmd = 'bedpeToBed12.py -i ' + bedpe + ' -d ' + str(dist) + ' > ' + bed12
    os.popen(cmd)
    
    bam = bed12 + ".bam"
    cmd = 'bedToBam -i ' + bed12 + ' -bed12 -g ' + genome + ' > ' + bam
    os.popen(cmd)
    
    sorted = bam + ".sorted"
    cmd = 'samtools sort ' + bam + ' ' + sorted
    os.popen(cmd)
    
    cmd = 'samtools index ' + sorted + '.bam'
    os.popen(cmd)
    
    out_idx = out + ".bai"
    shutil.move(sorted + ".bam", out)
    shutil.move(sorted + ".bam.bai", out_idx)
    
    os.remove(bed12)
    os.remove(bam)
    
def main():
	usage = """%prog -i <file> -n <name> -d <dist>

bedpeToBam version 1.0
Author: Aaron Quinlan & Ira Hall	
Description: converts BEDPE to sorted/index BAM format for viewing in IGV or the UCSC browser.
Last Modified: July 20, 2010 

	"""
	parser = OptionParser(usage)

	parser.add_option("-i", "--bedpe", dest="bedpe",
		help="BEDPE input file", metavar="FILE")

	parser.add_option("-d", "--maxdist", dest="dist", default = 100000, type="int",
		help="The minimum distance for drawing intrachromosomal features as if they are interchromosomal (i.e., without a line spanning the two footprints). Default is 1Mb.",
		metavar="INT")
		
	parser.add_option("-g", "--genome", dest="genome", type="string",
		help="The BEDTools genome file",
		metavar="STRING")
		
	parser.add_option("-o", "--out", dest="out", type="string",
		help="The output file name",
		metavar="STRING")

	(opts, args) = parser.parse_args()

	if opts.bedpe is None:
		parser.print_help()
		print
	else:
		processBEDPE(opts.bedpe, opts.dist, opts.out, opts.genome)

if __name__ == "__main__":
	main()



