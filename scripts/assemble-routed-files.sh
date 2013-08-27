if [ -z $3 ]
then
	echo "usage:$0 <file w/ list of routed files> <hydra config file> <number of threads>"
	exit
fi

CONFIG=$1
ROUTED_FILES=$2
PROCS=$3
COUNTER=0
INDEX=0

function add_next_assem {
    # if still jobs to do then add one
    if [[ $COUNTER -lt PROCS ]]
    then
	COUNTER=$(($COUNTER+1))
		hydra-assembler -config $CONFIG -routed ${assems_todo[$INDEX]} -maxMappings 300000 & 
		INDEX=$(($INDEX+1))
	else
		wait
		COUNTER=$(($COUNTER-$PROCS))
    fi
}

function add_final_assem {
    hydra-assembler -config $CONFIG -routed ${assems_todo[$INDEX]} -maxMappings 300000
    wait
}

assems_todo=($(cat $ROUTED_FILES)) # places output into an array
max_index=${#assems_todo[*]}-1
while [[ $INDEX -lt $max_index ]]
do
    add_next_assem
done
if [[ $INDEX -eq $max_index ]]
    then
    add_final_assem
fi
