#!/usr/bin/env python
import sys
import os
import pysam
from math import sqrt
from optparse import OptionParser


def get_median(l):
    """Calculate the median value of a list"""
    
    sorts = sorted(l)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]


def get_mad(l):
    """Calculate the median absolute deivation of a list"""
    
    median = get_median(l)
    # compute the MAD
    diff_from_median = []
    for val in l:
        diff_from_median.append(abs(val - median))
    return get_median(diff_from_median)


def meanstdv(x):
    """
    Calculate mean and standard deviation of data x[]:
        mean = {\sum_i x_i \over n}
        std = sqrt(\sum_i (x_i - mean)^2 \over n-1)
    """
    
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    return mean, std


def parse_config_stub(config_stub):
	"""Yield each sample and filename from stub config file by bam size."""
	stub_list = []
	for line in open(config_stub):
		fields = line.strip().split('\t')
		if len(fields) != 2:
			sys.exit("config stub file should only have 2 fields")
		sample, file = fields[0], fields[1]
		file_size = os.path.getsize(file)
		(sample, file, file_size) = sample, file, file_size
		stub_list.append((sample, file, file_size))
	"""Sort the stub list of sets by byte largest value"""
	stub_list.sort(key=lambda x: x[2], reverse=True)
	for item in stub_list:
		yield item


def main():
    usage = """%prog -i <config_stub>
    """
    parser = OptionParser(usage)

    parser.add_option("-i", dest="config_stub",
        help="Basic input sample file (sample_id[TAB]sample_file_path)",
        metavar="FILE")
    parser.add_option("-s", dest="sample_size",
        help="How many pairs to sample (def. 100000)",
        metavar="INT", type="int", default=100000)
    parser.add_option("-n", dest="num_variance_units",
        help="The num. of units of variation that should be allowed (def. 20)",
        metavar="INT", type="int", default=20)


    # Grab the command line options
    (opts, args) = parser.parse_args()

    # Make sure we have wha we need.  If so, on to the goodness.
    if opts.config_stub is None:
        parser.print_help()
        print
    else:
        for (sample, file, file_size) in parse_config_stub(opts.config_stub):

            bamfile = pysam.Samfile(file, "rb")
            isizes = []
            # collect representative insert sizes for this BAM file.
            for aln in bamfile:
                # only consider F/R pairs - ignore its R/F buddy in the file
                if not aln.is_reverse and aln.mate_is_reverse:
                    isizes.append(aln.tlen)
                    if len(isizes) >= opts.sample_size:
                        break
            bamfile.close()

            med = get_median(isizes)
            mad = get_mad(isizes)
            print '\t'.join([str(s) for s in \
                            [sample, file, med, mad, opts.num_variance_units]])

if __name__ == "__main__":
    main()



