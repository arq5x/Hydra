#Hydra-Multi - an SV discovery tool that incorporates hundreds of samples
=======================================================================

##Overview

Hydra-Multi is a paired-end read structural variant discovery tool that is capable of integrating signals from hundreds of samples.

##Installation

###Requirements
=========================
1. [samtools](http://samtools.sourceforge.net/)
2. [pysam](https://code.google.com/p/pysam/)
3. `set ulimit -f 16384`

The ulimit determines the number of open file handles on a system.  
This number must be larger than 4*number of possible chromosome-chromosome combinations found in the desired reference genome.  
For the human reference (hg19 b37), 16384 is the recommended ulimit.

####Installaing
	git clone https://github.com/arq5x/Hydra
	cd Hydra
	make 
	chmod +x Hydra/scripts/*
	sudo cp Hydra/scripts/* /usr/local/bin
	sudo cp Hydra/bin/* /usr/local/bin

####Testing Hydra-Multi
	chmod +x hydra-multi.sh
	./hydra-multi.sh test
	
##Running Hydra-Multi
==========================
A wrapper script (hydra-multi.sh) can be used to automatiically run Hydra-Multi or each step may be performed manually. Both the automatic and manual executions require a stub file to create a config file.  

0. Generate a stub file.
==========================
Start with a simple config file "stub" such as the one below:

    $ cat config.stub.txt
    sample1	/full/path/to/file/sample1.pos.bam
    sample2	/full/path/to/file/sample2.pos.bam
    sample3	/full/path/to/file/sample3.pos.bam

###Automatic Execution
==========================

hydra-multi.sh can then then be used to execute subsequent steps:

	./hydra-multi.sh run config.stub.txt

To obtain a parameter list for running hydra-multi:

	./hydra-multi.sh run -h

	
###Manual Execution 
1. Generate a config file.
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


2. Extract discordant alignments.
=================================
Once you have created a configuration file for HydraMulti, you need to run the
`extract_discordants.py` script to, you guessed it, extract the discordant 
alignments from your BAM files into BEDPE format for HydaMulti.

NOTE: the `extract_discordants.py` script inspects the is_proper_pair bit (0x2)
in the SAM FLAG field to identify discordant alignments.  If you want to use 
other rules for discordancy, you will need to write a script to set the FLAG
according to your custom rules.

By default, `extract_discordants.py` requires both ends of a pair to be aligned,
have MAPQ >= 20, and requires both ends of a pair to have an edit distance of at
most 4.  Moreover, it filters out any alignments marked as duplicates.  One can
override these settings with the `--min_mapq`, `--max_edit`, and `--allow_dups`
options, respectively.

For each inout BAM file listed in your configuration file, 
`extract_discordants.py` will create a BEDPE file of the discordant alignments
in the the same directory.  For example, it will create a `sample1.pos.bam.bedpe` 
file for the `sample1.pos.bam` input file listed in the config file:

    python scripts/extract_discordants.py -i config.hydra.txt


3. Run HydraRouter
=================================

    $ hydra-router -config config.hydra.txt -routedList routed-files.txt


4. Assemble SV breakpoint clusters
==================================

    $ sh scripts/assemble-routed-files.sh routed-files-test.txt config.hydra.txt


5. Combine the individual SV assembly files into a single file.
===============================================================

    $ sh scripts/combine-assembled-files.sh /full/path/to/assembled/files/ all.assembled


6. Finalize the SV breakpoint predictions.
===============================================================

    $ scripts/forceOneClusterPerPairMem.py -i all.assembled -o all.sv-calls
