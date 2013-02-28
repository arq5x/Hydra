if [ -z $2 ]
then
	echo "usage:$0 <<hydra config file> <number of processors>"
	exit
fi

CONFIG=$1
PROCS=$2

COUNTER=0
for bam in `cut -f 2 $CONFIG`; 
do 
	python scripts/extract_discordants.py -i $bam &
	COUNTER=$((COUNTER + 1))
	if (( $COUNTER % $PROCS == 0 )); then wait; fi # Limit to $PROCS concurrent subshells.
done