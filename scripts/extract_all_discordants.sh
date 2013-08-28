#!/bin/bash
set -e

if [ -z $2 ]
then
	echo "usage:$0 <<hydra config file> <number of processes>"
	exit
fi

CONFIG=$1
PROCS=$2
INDEX=0

function poll {
    while [[ $curr_jobs -ge PROCS ]]
    do
	curr_jobs=$(jobs -p | wc -l)
	sleep 2
	wait $!
    done
}


function add_next_bam {
    # if still jobs to do then add one
    curr_jobs=$(jobs -p | wc -l)
    if [[ $curr_jobs -lt PROCS ]]
    then
	extract_discordants ${bams_todo[$INDEX]} & 
	INDEX=$(($INDEX+1))
    else
	poll
    fi

}


function extract_discordants {
	samtools view -bF 0x040E -f 0x001 $1 > $1.disc.tmp.bam
	samtools sort -m 1000000000 -n $1.disc.tmp.bam $1.disc.tmp.bam.qrysort
	extract_discordants.py -i $1.disc.tmp.bam.qrysort.bam $2
}

dataset_names=($(cut -f 1 $CONFIG))
bams_todo=($(cut -f 2 $CONFIG)) # places output into an array
max_index=${#bams_todo[*]}-1
while [[ $INDEX -le $max_index ]]
do
    add_next_bam
done
wait



