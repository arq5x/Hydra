Hydra-Multi - an SV discovery tool that incorporates hundreds of samples
=======================================================================

#Overview

Hydra-Multi is a paired-end read structural variant discovery tool that is capable of integrating signals from hundreds to thousands of samples.

#Installation
Below are the requirements and instructions for installation of Hydra-Multi. 

##Requirements
1. [samtools](http://samtools.sourceforge.net/)
2. [pysam](https://code.google.com/p/pysam/)
3. `set ulimit -f 16384`

The ulimit determines the number of open file handles on a system.  
This number must be larger than 4*number of possible chromosome-chromosome combinations in the respective reference.  
For the human reference (hg19 b37), 16384 is the recommended ulimit.

###Installaing
	git clone https://github.com/arq5x/Hydra
	cd Hydra
	make 
	chmod +x scripts/*
	sudo cp scripts/* /usr/local/bin
	sudo cp bin/* /usr/local/bin

###Testing Hydra-Multi
	chmod +x hydra-multi.sh
	./hydra-multi.sh test
	
#Running Hydra-Multi
A wrapper script (hydra-multi.sh) can be used to automatiically run Hydra-Multi or each step may be performed manually. Both the automatic and manual executions begin by creating a stub file. 

###0. Generate a stub file.
==========================
Start with a simple config file "stub" such as the one below:

    $ cat config.stub.txt
    sample1	/full/path/to/file/sample1.pos.bam
    sample2	/full/path/to/file/sample2.pos.bam
    sample3	/full/path/to/file/sample3.pos.bam

##Automatic Execution
hydra-multi.sh can then then be used to execute subsequent steps:

	./hydra-multi.sh run config.stub.txt

To obtain a parameter list for using the run script:

	$./hydra-multi.sh run -h
	
	usage:   hydra-multi.sh run [options] <stub_file>
	
	positional arguments:
		stub file
			the stub file to create the configuration file, example on https://github.com/arq5x/Hydra
	options:
		-t INT	Number of threads to use. [Default: 2]
		-p INT	The punt parameter for assembly, the maximum read depth allowed. [Default: 10]
		-o STR	The stub for the output file names

	
##Manual Execution 

###1. Generate a config file.
==========================

HydraMulti needs a configuration file documenting the sample/libraries and the
paths to their respective BAM files that will be input to SV discovery process.

The `make_hydra_config.py` script will inspect the alignments in each sample's
BAM file to automatically create a complete config file documenting the
statistics of the fragment library:

    python scripts/make_hydra_config.py -i config.stub.txt
    sample1	/full/path/to/file/sample1.pos.bam	374.23	12	3
    sample2	/full/path/to/file/sample2.pos.bam	398.19	20	3
    sample3	/full/path/to/file/sample3.pos.bam	401.78	23	3
	
Just redirect the output to a new, complete config file and you should be
ready to go:

    python scripts/make_hydra_config.py -i config.stub.txt > config.hydra.txt


###2. Extract discordant alignments.
=================================
Once you have created a configuration file for Hydra-Multi, you need to run the
`extract_discordants.py` script to, you guessed it, extract the discordant 
alignments from your BAM files into BEDPE format for HydaMulti.

For each inout BAM file listed in your configuration file, 
`extract_discordants.py` will create a BEDPE file of the discordant alignments
in the the same directory.  For example, it will create a `sample1.pos.bam.bedpe` 
file for the `sample1.pos.bam` input file listed in the config file:

    python scripts/extract_discordants.py -c config.hydra.txt -d <sample_name>


###3. Run HydraRouter
=================================
This routes all of the alignments on with the same chromosome/orientation set to the same file for assembly.

    $ hydra-router -config config.hydra.txt -routedList routed-files.txt


###4. Assemble SV breakpoint clusters
==================================
Assembly of each chromosome/orientation set.

    $ sh scripts/assemble-routed-files.sh routed-files-test.txt config.hydra.txt 1


###5. Combine the individual SV assembly files into a single file.
===============================================================
Combine all of the chromosome/orientation sets back into one file.

    $ sh scripts/combine-assembled-files.sh /full/path/to/assembled/files/ all.assembled


###6. Finalize the SV breakpoint predictions.
===============================================================

    $ scripts/forceOneClusterPerPairMem.py -i all.assembled -o all.sv-calls


###7. Determine presence of the SV breakpoint predictions in samples.
=======================================================================

    $ scripts/frequency.py -f all.sv-calls.final -d all.sv-calls.detail > all.sv-calls.freq
    
    
###8. Change footprint intervals into breakpoint intervals.
===============================================================

    $ scripts/hydraToBreakpoint -i all.sv-calls.freq >  all.sv-calls.bkpts
    
