#!/usr/bin/env python
import collections
import sys
import os
import pysam
import string
import subprocess
import shlex
import multiprocessing as mp
from optparse import OptionParser

class BEDPE (object):

    def __init__(self, bam1, bam2, bamfilehandle, datasetname):
        
        if ((bam1.tid == bam2.tid and bam2.pos < bam1.pos) or \
           (bam1.tid != bam2.tid and bam2.tid < bam1.tid)):
            bam2, bam1 = bam1, bam2
                
        self.c1        = bamfilehandle.getrname(bam1.tid)
        self.s1        = bam1.pos
        self.e1        = bam1.aend
        self.c2        = bamfilehandle.getrname(bam2.tid)
        self.s2        = bam2.pos
        self.e2        = bam2.aend
        self.name      = datasetname + "." + bam1.qname
        
        if bam1.is_read1:
            self.score     = 1
        else:
            self.score     = 2
            
        self.o1        = "-" if bam1.is_reverse else "+"
        self.o2        = "-" if bam2.is_reverse else "+"
        self.edit1     = self.get_edit_distance(bam1)
        self.edit2     = self.get_edit_distance(bam2)
        self.mappings1 = 1   # used to prioritize deduping.
        self.mappings2 = 1   # used to prioritize deduping.
        self.mapq1     = int(bam1.mapq)
        self.mapq2     = int(bam2.mapq)
        self.cigar1    = bam1.cigar
        self.cigar2    = bam2.cigar
        self.unmap1    = bam1.is_unmapped
        self.unmap2    = bam2.is_unmapped
        

    def get_edit_distance(self, bam):
        nm_tag = [tag for tag in bam.tags if tag[0] == "NM"]
        if len(nm_tag) >= 1:
            return nm_tag[0][1]
        elif len(nm_tag) == 0:
            return 0

    def __repr__(self):
        return '\t'.join([str(s) for s in [self.c1, self.s1, self.e1,
                                          self.c2, self.s2, self.e2, 
                                          self.name, self.score,
                                          self.o1, self.o2,
                                          self.edit1, self.edit2,
                                          self.mappings1, self.mappings2,
                                          self.mapq1, self.mapq2]])


def parse_config(config, dataset):
    """Yield each sample and filename and stats from  config file."""
    
    for line in open(config):
        fields = line.strip().split('\t')
        if len(fields) != 5:
            sys.exit("config file should only have 5 fields")
        if (fields[0] == dataset):
            sample = fields[0]
            bamFile = fields[1]
            centralStat = float(fields[2])
            deviation = float(fields[3])
            numDevs = float(fields[4])
        else:
            continue
    if sample == None:
        sys.exit("config file does not contain dataset: " + dataset + "\n")
    return sample, bamFile, centralStat, deviation, numDevs
 

def make_discordant_bam(bam, allow_dups, sampleName, expected, dev, numDev):
    bamfile = pysam.Samfile(bam, "rb")
    tmp_bam_file_name = bam + ".disc.tmp.bam"
    tmp_bamfile = pysam.Samfile(tmp_bam_file_name, "wb", template = bamfile)
    outDist = dev*numDev
    aln = None
    for aln in bamfile:
        #toss unmapped/orphans: samtools view -uF 0x0004 bam 
       	#| samtools view -uF 0x0008 - | samtools view -uf 0x0001 
        if (not aln.is_unmapped and not aln.mate_is_unmapped and aln.is_paired):
			#is duplicate?
            if (aln.is_duplicate and not allow_dups):
                continue
            #interchrom auto-disc
            if (aln.tid != aln.rnext):
                tmp_bamfile.write(aln)
                #-/-
            elif (aln.is_reverse and aln.mate_is_reverse):
                tmp_bamfile.write(aln)
                #+/+
            elif (not aln.is_reverse and not aln.mate_is_reverse):
                tmp_bamfile.write(aln)
                #is insert/deletion? : ((on same chrome) and ((r1+/r2-) or (r2+/r2-)))?
            elif ((aln.tid == aln.rnext) and (not aln.is_reverse and aln.mate_is_reverse and (int(aln.pos)<=int(aln.pnext))) or \
			(aln.is_reverse and not aln.mate_is_reverse and (int(aln.pos)>=int(aln.pnext)))):
				#is the alignment a delete?
                if (abs(aln.tlen) > expected + outDist):
                    tmp_bamfile.write(aln)
				#is the alignment an insert?
                #elif (abs(aln.tlen) < expected - outDist):
                #    tmp_bamfile.write(aln)
                else:
                        #in concordant range
                    continue
            #this discordant read supports a non-deletion/insert type event
            else:
                tmp_bamfile.write(aln)
    tmp_bamfile.close()
    bamfile.close()
    return tmp_bam_file_name

def query_sort_discordant(bam_filename, mem):
    query_sort_filename = bam_filename + ".qrysort"
    sortCmd = "samtools sort -n -m " + str(mem) + " " + bam_filename + " " + query_sort_filename
    subprocess.call(shlex.split(sortCmd))
    os.remove(bam_filename + ".disc.tmp.bam)
    return query_sort_filename + ".bam"

def make_discordant_bedpe(discordant_bam_filename, 
	                          min_mapq, dataset_name):
	
	orig_bam_filename = discordant_bam_filename
	idx = string.find(orig_bam_filename, ".disc.tmp.bam.qrysort.bam")
	bedpe_outfile_name = orig_bam_filename[0:idx] + ".bedpe"
	bedpe_outfile = open(bedpe_outfile_name, 'w')
	bam1 = None
	bam2 = None
	bamfile = pysam.Samfile(discordant_bam_filename, "rb")    
	for aln in bamfile:
		bam1 = aln
		try:
			bam2 = bamfile.next()
		except StopIteration:
			break
		if (bam1.qname != bam2.qname):
			while (bam1.qname != bam2.qname):
				bam1 = bam2
				try:
					bam2 = bamfile.next()
				except StopIteration:
					break
		elif bam1.is_paired and bam2.is_paired and \
		not bam1.is_unmapped and not bam2.is_unmapped :
			pair = BEDPE(bam1, bam2, bamfile, dataset_name)
			
			if int(pair.mapq1) < int(min_mapq) or int(pair.mapq2) < int(min_mapq):
				continue
	        
			bedpe_outfile.write(str(pair) + '\n')
	
	bedpe_outfile.close()
    os.remove(discordant_bam_filename)
	return bedpe_outfile


def main():
    usage = """%prog -c <config file> -d <dataset>
    """
    parser = OptionParser(usage)
    parser.add_option("-c", dest="config_file",
		help="Input config file (Required)",
		metavar="FILE")
    parser.add_option("-d", dest="dataset",
		help="Dataset name (Required)",
		metavar="STRING")
    parser.add_option("--min_mapq", dest="min_mapq",
        help="Minimum MAPQ required on both ends of a pair (def. 0)",
        metavar="INT", default=0)  
    parser.add_option('--allow_dups', dest="allow_dups",
        help="Allow duplicate pairs (def. False)",
        action="store_true", default=False)
    parser.add_option('--mem', dest="mem",
	    help="Memory used for sorting (def. 2000000000)",
		metavar="INT", default=2000000000)
    # Grab the command line options
    (opts, args) = parser.parse_args()

    # Make sure we have wha we need.  If so, on to the goodness.
    if (opts.config_file or opts.dataset) is None:
        parser.print_help()
        print
    else:
        (sample, bamFile, expec, dev, numDev)  = parse_config(opts.config_file, opts.dataset)
        #pull out our reads to make a tmp bam
       	discordant_bam = make_discordant_bam(bamFile,
		opts.allow_dups,
		sample,
		expec,
		dev,
		numDev)
        discordant_bam_query_sort = query_sort_discordant(discordant_bam, opts.mem)
		
        discordant_bedpe = make_discordant_bedpe(discordant_bam_query_sort,
		opts.min_mapq,
		sample)


if __name__ == "__main__":
    main()



