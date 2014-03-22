#!/bin/bash
usage()
{
	echo "
	usage: hydra-multi.sh <command> [options]
	command:
	test	perform a test run of Hydra-Multi with 3 sample datasets
	run	perform a run of Hydra-Multi
	options: -h	show this help
"
}

function test() {
	PUNT=10
	THREADS=2
	echo "Downloading 3 sample files from 1000 Genomes (~1.5Gb total)...\c"
	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00689/alignment/HG00689.chrom11.ILLUMINA.bwa.CHS.low_coverage.20120522.bam
	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG01615/alignment/HG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam
	echo "done"
	
	
	echo "creating a basic configuration file from the downloaded 1000G files...\c"
	echo -e "HG00096\tHG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
	HG00689\tHG00689.chrom11.ILLUMINA.bwa.CHS.low_coverage.20120522.bam
	HG01615\tHG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam" > config.stub.txt
	echo "done"
	
	
	echo "creating a complete configuration file by sampling BAM to create library stats...\c"
	python scripts/make_hydra_config.py -i config.stub.txt > config.hydra.txt
	echo "done"
	
	
	echo "extracting discordant alignments from BAM files...\c"
	sh scripts/extract_all_discordants.sh config.hydra.txt $THREADS
	echo "done"
	
	
	echo "running hydra router on the discordant alignments...\c"
	hydra-router -config config.hydra.txt -routedList routed-files.txt
	echo "done"
	
	
	echo "running hydra-assembler on the routed files of discordant alignments...\c"
	sh scripts/assemble-routed-files.sh config.hydra.txt routed-files.txt $THREADS $PUNT
	echo "done"
	
	
	echo "re-combining the individual assembled files...\c"
	sh scripts/combine-assembled-files.sh   ./   all.1000G.assembled
	echo "done"
	
	
	echo "finalizing SV breakpoint calls...\c"
	python scripts/finalizeBreakpoints.py -i all.1000G.assembled -o all.1000G.sv
	echo "done"
}


#function run() {
#	function run_usage() {
#		"echo how to use test
#		"
#	}
#	if test -z "$1"; then
#		run_usage
#		exit 1
#    	fi
#	echo "Downloading 3 sample files from 1000 Genomes (~1.5Gb total)...\c"
#	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
#	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00689/alignment/HG00689.chrom11.ILLUMINA.bwa.CHS.low_coverage.20120522.bam
#	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG01615/alignment/HG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam
#	echo "done"
#	
#	
#	echo "creating a basic configuration file from the downloaded 1000G files...\c"
#	echo -e "HG00096\tHG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
#	HG00689\tHG00689.chrom11.ILLUMINA.bwa.CHS.low_coverage.20120522.bam
#	HG01615\tHG01615.chrom11.ILLUMINA.bwa.IBS.low_coverage.20120522.bam" > config.stub.txt
#	echo "done"
#	
#	
#	echo "creating a complete configuration file by sampling BAM to create library stats...\c"
#	python scripts/make_hydra_config.py -i config.stub.txt > config.hydra.txt
#	echo "done"
#	
#	
#	echo "extracting discordant alignments from BAM files using 2 threads...\c"
#	sh scripts/extract_all_discordants.sh config.hydra.txt 2
#	echo "done"
#	
#	
#	echo "running hydra router on the discordant alignments...\c"
#	hydra-router -config config.hydra.txt -routedList routed-files.txt
#	echo "done"
#	
#	
#	echo "running hydra-assembler on the routed files of discordant alignments using 2 threads and punting at read depth of 10...\c"
#	sh scripts/assemble-routed-files.sh config.hydra.txt routed-files.txt $THREADS $PUNT
#	echo "done"
#	
#	
#	echo "re-combining the individual assembled files...\c"
#	sh scripts/combine-assembled-files.sh   ./   all.1000G.assembled
#	echo "done"
#	
#	
#	echo "finalizing SV breakpoint calls...\c"
#	python scripts/finalizeBreakpoints.py -i all.1000G.assembled -o all.1000G.sv
#	echo "done"
#}

if test -z "$1"
then
    usage
    exit 1
fi

while getopts "h" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

case "$1" in 
    'test')
	test
	;;
    *)
	usage
	echo -e "Error: command \"$1\" not recognized\n"
	exit 1
esac

