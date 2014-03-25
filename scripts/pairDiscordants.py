#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser

def pairReads(inFile, numMappings, order, dist, minSpan, minConcRange, maxConcRange, mode, anchorThresh, multiThresh, editSlop, writeOrphans):
	lastID = 0
	currentID = 0
	xval = []
	yval = []
	read1 = []
	read2 = []
	mappings1 = 0
	mappings2 = 0
	badID = 0
	if inFile == "stdin":
		data = sys.stdin
	else:
		data = open(inFile, 'r')
	while 1:
		line = data.readline()
		if not line:
			pairs = makePairs(read1, read2, numMappings, order, dist, minSpan, minConcRange, maxConcRange, mode, mappings1, mappings2)
			if pairs != "BAD":
				if mode == "hydra":
					printHydraMappings(pairs, editDistance, editSlop, False)
				elif mode == "discordant":
					printDiscordantMappings(pairs, editDistance, editSlop)
				elif mode == "detail":
					printDetailMappings(pairs, editDistance, editSlop)
			break
			
		xval = line.strip().split('\t')
		idSplit = xval[3].split('/')
		yval = [str(xval[0]), str(xval[1]), str(xval[2]), str(idSplit[0]), str(idSplit[1]), str(xval[4]), str(xval[5])]
		currentID = yval[3] 
		read = int(yval[4])
		if currentID == badID:
			continue
		if ((mappings1 == 1) and (mappings2 > anchorThresh)) or ((mappings2 == 1) and (mappings1 > anchorThresh)) or ((mappings1 > 1) and (mappings2 > 1) and ((mappings1 * mappings2) > multiThresh)):
			badID = lastID
		if (currentID == lastID) or (lastID == 0):
			lastID = currentID
			if read == 1:
				read1.append(yval)
				mappings1 = mappings1 + 1
			elif read == 2: 
				read2.append(yval)	
				mappings2 = mappings2 + 1
		elif currentID != lastID:
			if ((mappings1 == 1) and (mappings2 > anchorThresh)) or ((mappings2 == 1) and (mappings1 > anchorThresh)) or ((mappings1 > 1) and (mappings2 > 1) and ((mappings1 * mappings2) > multiThresh)):
				#print "THIS IS A REPEAT"
				if mode == "detail":
					print "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + lastID + "\t" + "NA" + "\t" + \
						"NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "repeat"				
			elif (mappings1 > 0) and (mappings2 > 0):
				pairs = makePairs(read1, read2, numMappings, order, dist, minSpan, minConcRange, maxConcRange, mode, mappings1, mappings2)
				if pairs != "BAD":
					editDistance = pairs[0][15]
					#print str(edit)
					if mode == "hydra":
						printHydraMappings(pairs, editDistance, editSlop, False)
					elif mode == "discordant":
						printDiscordantMappings(pairs, editDistance, editSlop)
					elif mode == "detail":
						printDetailMappings(pairs, editDistance, editSlop)
			else:
				if mode == "detail":
					print "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + lastID + "\t" + "NA" + "\t" + \
					"NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "orphan"
				elif ((mode == "hydra") and (writeOrphans == True)):
					printHydraMappings(pairs, 0, editSlop, True)
							
			xval = line.strip().split('\t')
			idSplit = xval[3].split('/')
			yval = [str(xval[0]), str(xval[1]), str(xval[2]), str(idSplit[0]), str(idSplit[1]), str(xval[4]), str(xval[5])]
			lastID = yval[3]
			read1 = []
			read2 = []
			if read == 1:
				read1.append(yval)
				mappings1 = 1
				mappings2 = 0
			elif read == 2: 
				read2.append(yval)	
				mappings1 = 0
				mappings2 = 1


def makePairs(read1, read2, numMappings, order, dist, minSpan, minConcRange, maxConcRange, mode, mappings1, mappings2):
	output = []
	badout = []
	new = []
	stop = 0
	count1 = 0
	count2 = 0
	for line in read1:
		if stop == 1:
			break
		count1 = count1 + 1
		count2 = 0
		for row in read2:
			count2 = count2 + 1
			chrom1 = line[0]
			pos1 = int(line[1])
			pos2 = int(line[2])
			strand1 = line[6]
			chrom2 = row[0]
			pos3 = int(row[1])
			pos4 = int(row[2])
			strand2 = row[6]
			readLength1 = pos2 - pos1
			readLength2 = pos4 - pos3
			nonOverlap = max(pos1, pos3) - min(pos2, pos4)	
			fragSize = nonOverlap + (readLength1 + readLength2)
			if (chrom1 != chrom2) or (fragSize > 1000000):
				fragSize = "distant" 
			if ((chrom1 == chrom2) and (strand1 == "+") and (strand2 == "-") and ((pos4 - pos1) >= minConcRange) and ((pos4 - pos1) < maxConcRange)) or \
				((chrom1 == chrom2) and (strand1 == "-") and (strand2 == "+") and ((pos2 - pos3) >= minConcRange) and ((pos2 - pos3) < maxConcRange)):
				newBad = [row[0], row[1], row[2], line[0], line[1], line[2], line[3], row[4], row[6], line[6], row[5], line[5], fragSize, mappings1, mappings2,	 int(line[5]) + int(row[5]), count1 + count2, "concordant"]
				stop = 1 
				break
			elif (chrom1 == chrom2) and (nonOverlap <= minSpan):
				newBad = [row[0], row[1], row[2], line[0], line[1], line[2], line[3], row[4], row[6], line[6], row[5], line[5], fragSize, mappings1, mappings2,	 int(line[5]) + int(row[5]), count1 + count2, "badSpan"]
				stop = 1
				break
			elif (chrom1 == chrom2) and (pos3 < pos1 and (fragSize != "distant")):
				new = [row[0], row[1], row[2], line[0], line[1], line[2], line[3], row[4], row[6], line[6], row[5], line[5], fragSize, mappings1, mappings2,  int(line[5]) + int(row[5]), count1 + count2, "discordant"]
				output.append(new)
			else:
				new = [line[0], line[1], line[2], row[0], row[1], row[2], row[3], line[4], line[6], row[6], line[5], row[5], fragSize, mappings1, mappings2, int(line[5]) + int(row[5]), count1 + count2, "discordant"]
				output.append(new)	
	# first sort by combined counts or size:
	if order == "genome":
		output.sort(cmp=lambda x,y: cmp(x[16],y[16]))
	elif order == "size":
		output.sort(cmp=lambda x,y: cmp(x[12],y[12]))
	# then by combined edit distance
	output.sort(cmp=lambda x,y: cmp(x[15],y[15]))
	if (stop == 1):
		if (mode == "discordant") or (mode == "hydra"):
			final = 'BAD'
		elif mode == "detail":
			badout.append(newBad)
			final = badout
	elif len(output) > numMappings:
		final = output[0:numMappings]
	else:
		final = output
		#print str(output[0][15])
	return(final)


# make the print statement conditional on the combined edit distance and range (e.g., within 2); pass the edit distance as an output right after the last sort; use the edit distance as input to print; 
def printHydraMappings(mappings, editDistance, editSlop, isOrphan):
	try:
		if isOrphan == False:
			for line in mappings:
				if (line[15] - editDistance) <= editSlop:
					if (str(line[0]) > str(line[3])) or ((str(line[0]) == str(line[3])) and (str(line[1]) > str(line[4]))):
						print str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) + "\t" + str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) \
						+ "\t" + str(line[6]) + "\t" + "2" + "\t" + str(line[9]) + "\t" + str(line[8]) + "\t" + str(line[11]) + "\t" + str(line[10])  + "\t" + str(line[14])  + "\t" + str(line[13]) 
					else:
						print str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) \
						+ "\t" + str(line[6]) + "\t" + "1" + "\t" + str(line[8]) + "\t" + str(line[9]) + "\t" + str(line[10]) + "\t" + str(line[11])  + "\t" + str(line[13])  + "\t" + str(line[14]) 
		else:
		    for line in mappings:
			    print str(line)
			
	except IOError, e:
		if e.errno != 32:
			raise	


def printDiscordantMappings(mappings, editDistance, editSlop):
	try:
		for line in mappings:
			if (line[15] - editDistance) <= editSlop:
				if (str(line[0]) > str(line[3])) or ((str(line[0]) == str(line[3])) and (str(line[1]) > str(line[4]))):
					print str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) + "\t" + str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) \
					+ "\t" + str(line[6]) + "\t" + "2" + "\t" + str(line[9]) + "\t" + str(line[8]) + "\t" + str(line[11]) + "\t" + str(line[10]) \
					+ "\t" + str(line[12])	+ "\t" + str(line[14])	+ "\t" + str(line[13])	+ "\t" + str(line[15])
				else:
					print str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) \
					+ "\t" + str(line[6]) + "\t" + "1" + "\t" + str(line[8]) + "\t" + str(line[9]) + "\t" + str(line[10]) + "\t" + str(line[11]) \
					+ "\t" + str(line[12])	+ "\t" + str(line[13])	+ "\t" + str(line[14])	+ "\t" + str(line[15])
	except IOError, e:
		if e.errno != 32:
			raise


def printDetailMappings(mappings, editDistance, editSlop):
	try:
		for line in mappings:
			if (line[15] - editDistance) <= editSlop:
				if (str(line[0]) > str(line[3])) or ((str(line[0]) == str(line[3])) and (str(line[1]) > str(line[4]))):
					print str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) + "\t" + str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) \
					+ "\t" + str(line[6]) + "\t" + "2" + "\t" + str(line[9]) + "\t" + str(line[8]) + "\t" + str(line[11]) + "\t" + str(line[10]) \
					+ "\t" + str(line[12])	+ "\t" + str(line[14])	+ "\t" + str(line[13])	+ "\t" + str(line[15])	+ "\t" + str(line[16])	+ "\t" + str(line[17])
				else:
					print str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) \
					+ "\t" + str(line[6]) + "\t" + "1" + "\t" + str(line[8]) + "\t" + str(line[9]) + "\t" + str(line[10]) + "\t" + str(line[11]) \
					+ "\t" + str(line[12])	+ "\t" + str(line[13])	+ "\t" + str(line[14])	+ "\t" + str(line[15])	+ "\t" + str(line[16])	+ "\t" + str(line[17])
	except IOError, e:
		if e.errno != 32:
			raise


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg		


def main(): 
	usage = """%prog -i <file> -n <numMappings, default=100> -o <order, default="size"> -d <distance, default=1000000> -l <readLength, default=42> -x <minSpan, default=100> -y <minConcLength, default=0> -z <maxConcLenth, default=600> -m <mode, default="discordant">

pairDiscordants version 1.0
Author: Ira Hall	
Description: Pairs and filters reads in BED format to make BEDPE files suitable for paired-end mapping with Hydra. 
Last Modified: July 20, 2010 
	
# INPUT: 
# chrom; start; end; readID/readNum; editDistance; strand
# chr10 14800035	14800077	113100100400208/1	0	+

# OUTPUT
# 1. chrom1
# 2. start1
# 3. end1
# 4. chrom2
# 5. start2
# 6. end2
# 7. read ID
# 8. to which end of the pair do fields 1,2,3,9,11 correspond? (1 or 2)
# 9. strand1 (+/-)
# 10. strand2
# 11. editDistance1
# 12. editDistance2		
# 13: mappings1
# 14: mappings2
# 15: fragSize (in discordant and detail mode only)
# 16: combined edit distance (in discordant and detail mode only)
# 17: summed count of row each mapping was present in block of original file (only in detail mode)
# 18: tag describing whether a read-pair was classified as discordant, concordant, badSpan, repeat or orphan (only in detail mode)
	"""
	parser = OptionParser(usage)
	
	parser.add_option("-i", "--file", dest="inFile", 
		help="A tab delimited BED file, or standard input (-i stdin).",
		metavar="FILE")
	
	parser.add_option("-n", "--numMappings", dest="numMappings", default=1000, type = "int",
		help="The number of mapping combinations that are reported. Default = 1000.",
		metavar="INT")
	
	parser.add_option("-o", "--order", dest="order", default="size",
		help="The order to choose mappings when there are more than specified in -n. Supports size between read1-read2 mappings (-o size) or genome order (-o genome) sorted in default fashion (chr10 before chr2). In either case mapping combinations with less edit distance have priority. Default = size.",
		metavar="STR")		
		
	parser.add_option("-d", "--dist", dest="dist", default=1000000, type = "int",
		help="The genomic distance above which a readpair mapping is considered 'distant'. This is used in determining the order that mappings are reported when using the 'size' option in -o. Default = 1000000",
		metavar="INT")

	parser.add_option("-x", "--minSpan", dest="minSpan", default=100, type = "int",
		help="The minimum genomic distance that a readpair must span in order to be reported, not counting the read lengths (i.e., inner span, not fragment size). Default = 100.",
		metavar="INT")

	parser.add_option("-y", "--minConcRange", dest="minConcRange", default=0, type = "int",
		help="minimum size to judge a readpair as 'concordant' with respect to the reference genome. Assumes fragment libraries not matepair libraries (+/- orientation is concordant). Default = 0.",
		metavar="INT")

	parser.add_option("-z", "--maxConcRange", dest="maxConcRange", default=600, type = "int",
		help="Maximum size to judge a read as 'concordant. Default = 600.",
		metavar="INT")

	parser.add_option("-m", "--mode", dest="mode", default="hydra",
		help="Report all mapping types (-m detail), just discordants (-m discordant), or directly to Hydra input format (-m hydra) without additional columns. In detail mode only one mapping is returned for concordant reads or those that fail the minSpan (-s1) threshold. These are tagged in an additional column as 'concordant' or 'badSpan'. Default = 'hydra'.",
		metavar="STR")

	parser.add_option("-a", "--anchorThresh", dest="anchorThresh", default=1000, type = "int",
		help="The number of mapping combinations allowed for readpairs with a single uniquely-mapped read; readpairs with more than this number of mappings on the non-unique end will not be reported (except in detail mode with NA's in all fields except the read ID); Default = 1000.",
		metavar="INT")
	
	parser.add_option("-r", "--multiThresh", dest="multiThresh", default=100, type = "int",
		help="The number of mapping combiantions allowed for readpairs without a single uniquely-mapped read; readpair mappings with more than this number of mappings will not be reported (except in detail mode with NA's in all fields except the read ID); Default = 100.",
		metavar="INT")

	parser.add_option("-e", "--editSlop", dest="editSlop", default=0, type = "int",
		help="The combined edit distance allowed, relative to best mapping for each pair; Default = 0.",
		metavar="INT")
		
	parser.add_option("-p", "--orphans", dest="writeOrphans", default=False, action="store_true",
		help="Whether or not to report \"orphaned\" pairs; Default = false.",
		metavar="INT")
			
	(opts, args) = parser.parse_args()

	if opts.inFile is None:
		parser.print_help()
		print
	else:
		pairReads(opts.inFile, opts.numMappings, opts.order, 
				  opts.dist, opts.minSpan, opts.minConcRange, 
				  opts.maxConcRange, opts.mode, opts.anchorThresh, 
				  opts.multiThresh, opts.editSlop, opts.writeOrphans)

#############################################
if __name__ == "__main__":
	sys.exit(main()) 