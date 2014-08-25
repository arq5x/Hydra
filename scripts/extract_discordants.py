#!/usr/bin/env python
import sys
import os
import pysam
import string
import subprocess
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


def parse_config(config):
    """Yield each sample and filename from  config file."""
    for line in open(config):
        fields = line.strip().split('\t')
        if len(fields) != 5:
            sys.exit("config stub file should only have 5 fields")
        
        (sample, file) = fields[0], fields[1]
        yield (sample, file)


# def make_discordant_bam(bam_filename):
# 
#     bamfile = pysam.Samfile(bam_filename, "rb")
# 
#     tmp_bam_file_name = bam_filename + ".disc.tmp.bam"
#     tmp_bamfile = pysam.Samfile(tmp_bam_file_name, "wb", 
#                                 template = bamfile)
# 
#     
#     # extract only the discordant BAM files.
#     for aln in bamfile:
#         if not aln.is_proper_pair:
#             tmp_bamfile.write(aln)
#     tmp_bamfile.close()
#     bamfile.close()
#     
#     return tmp_bam_file_name


# def query_sort_discordant(bam_filename):
#     query_sort_filename = bam_filename + ".qrysort"
#     pysam.sort("-m", "1000000000", "-n", 
#                bam_filename, query_sort_filename )
#     os.remove(bam_filename)
#     return query_sort_filename + ".bam"
# 

def make_discordant_bedpe(discordant_bam_filename, 
                          min_mapq, max_edit, allow_dups, dataset_name):
    
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
            
            if int(pair.edit1) > int(max_edit) or int(pair.edit2) > int(max_edit):
                continue
            
            if (bam1.is_unmapped and not allow_dups) or \
               (bam2.is_duplicate and not allow_dups):
                continue

            bedpe_outfile.write(str(pair) + '\n')
    
    bedpe_outfile.close()
    os.remove(discordant_bam_filename)
    return bedpe_outfile


def main():
    usage = """%prog -i <bam_file> -d <dataset_name>
    """
    parser = OptionParser(usage)

    parser.add_option("-i", dest="bam_file",
        help="Input BAM file",
        metavar="FILE")
        
    parser.add_option("-d", dest="dataset",
        help="Dataset name",
        metavar="FILE")


    
    parser.add_option("--min_mapq", dest="min_mapq",
        help="Minimum MAPQ required on both ends of a pair (def. 20)",
        metavar="INT", default=20)
        
    parser.add_option("--max_edit", dest="max_edit",
        help="Maximum edit distance allowed on both ends of a pair (def. 4)",
        metavar="INT", default=4)
        
    parser.add_option('--allow_dups', dest="allow_dups",
        help="Allow duplicate pairs (def. False)",
        action="store_true", default=False)


    # Grab the command line options
    (opts, args) = parser.parse_args()

    # Make sure we have wha we need.  If so, on to the goodness.
    if opts.bam_file is None:
        parser.print_help()
        print
    else:
        #discordant_bam = make_discordant_bam(opts.bam_file)
        # discordant_bam_query_sort = query_sort_discordant(discordant_bam)
        discordant_bedpe = make_discordant_bedpe(opts.bam_file,
                                                 opts.min_mapq,
                                                 opts.max_edit,
                                                 opts.allow_dups,
                                                 opts.dataset)


if __name__ == "__main__":
    main()



