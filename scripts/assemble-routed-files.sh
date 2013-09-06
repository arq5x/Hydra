#!/bin/bash
set -e

if [ -z $3 ]
then
	echo "usage:$0 <file w/ list of routed files> <hydra config file> <number of processes>"
	exit
fi

CONFIG=$1
ROUTED_FILES=$2
PROCS=$3
INDEX=0

function poll {
    while [[ $curr_jobs -ge PROCS ]]
    do
	curr_jobs=$(jobs -p | wc -l)
	sleep 1 &
	wait $!
    done
}

function add_next_assem {
    # if still jobs to do then add one
    curr_jobs=$(jobs -p | wc -l)
    if [[ $curr_jobs -lt PROCS ]]
    then
	hydra-assembler -config $CONFIG -routed ${assems_todo[$INDEX]} -maxMappings 1000*${#assems_todo[*]} & 
	INDEX=$(($INDEX+1))
    else
	poll
    fi
}


assems_todo=($(cat $ROUTED_FILES)) # places output into an array
max_index=${#assems_todo[*]}-1
while [[ $INDEX -le $max_index ]]
do
    add_next_assem
done
wait
