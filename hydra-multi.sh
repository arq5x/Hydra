#!/bin/bash
function usage()
{
	echo "
	usage: hydra-multi.sh <command> [options]
	command:
	test	perform a test of Hydra-Multi using 3 sample 1000 Genomes datasets
	run	execute Hydra-Multi
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
	python make_hydra_config.py -i config.stub.txt > config.hydra.txt
	echo "done"
	
	
	echo "extracting discordant alignments from BAM files...\c"
	sh extract_all_discordants.sh config.hydra.txt $THREADS
	echo "done"
	
	
	echo "running hydra router on the discordant alignments...\c"
	hydra-router -config config.hydra.txt -routedList routed-files.txt
	echo "done"
	
	
	echo "running hydra-assembler on the routed files of discordant alignments...\c"
	sh assemble-routed-files.sh config.hydra.txt routed-files.txt $THREADS $PUNT
	echo "done"
	
	
	echo "re-combining the individual assembled files...\c"
	sh combine-assembled-files.sh   ./   all.1000G.assembled
	echo "done"
	
	
	echo "finalizing SV breakpoint calls...\c"
	python finalizeBreakpoints.py -i all.1000G.assembled -o all.1000G.sv
	echo "done"
}


function run() {
	function run_usage() {
		echo "
	usage:   hydra-multi.sh run [options] <stub_file>
	
	positional arguments:
		stub file
			the stub file to create the configuration file, example on https://github.com/arq5x/Hydra
	options:
		-t INT	Number of threads to use. [Default: 2]
		-p INT	The punt parameter for assembly, the maximum read depth allowed. [Default: 10]
		-o STR	The stub for the output file names"
	}
	
	if [[ -z "$1" ]]; then
		run_usage
		exit 1
	fi
	THREADS=2
	PUNT=2
	OUT="hydra"
	STUB="${@:${OPTIND}:1}"
	while getopts ":t:p:o:" OPTION
	do
		case "${OPTION}" in
			h)
				run_usage
				exit 1
				;;
			t)
				THREADS="$OPTARG"
				;;
			o)
				OUT="$OPTARG"
				;;
		esac
	done
	
	echo "creating a complete configuration file by sampling BAM to create library stats...\c"
	python make_hydra_config.py -i $STUB > config.$OUT.txt
	echo "done"
	
	
	echo -e "extracting discordant alignments from BAM files using " $THREADS " threads...\c"
	sh extract_all_discordants.sh config.$OUT.txt $THREADS
	echo "done"
	
	
	echo "running hydra router on the discordant alignments...\c"
	hydra-router -config config.$OUT.txt -routedList routed-files.$OUT.txt
	echo "done"
	
	
	echo -e "running hydra-assembler on the routed files of discordant alignments using " $THREADS " threads and punting at read depth of 10...\c"
	sh assemble-routed-files.sh config.$OUT.txt routed-files.$OUT.txt $THREADS $PUNT
	echo "done"
	
	
	echo "re-combining the individual assembled files...\c"
	sh combine-assembled-files.sh   ./   all.$OUT.assembled
	echo "done"
	
	
	echo "finalizing SV breakpoint calls...\c"
	python finalizeBreakpoints.py -i all.$OUT.assembled -o all.$OUT.sv
	echo "done"
}

if [ -z $1 ]; then
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
			exit 1
			;;
	esac
done

case "$1" in 
	'run')
		run "${@:2}"
		;;
	'test')
		test
		;;
	*)
		usage
		echo -e "Error: command \"$1\" not recognized\n"
		exit 1
esac

