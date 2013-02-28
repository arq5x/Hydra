if [ -z $2 ]
then
	echo "usage:$0 <file w/ list of routed files> <hydra config file>"
	exit
fi

CONFIG=$1
ROUTED_FILES=$2


for routed in `cat $ROUTED_FILES`; 
    do hydra-assembler -config $CONFIG -routed $routed -maxMappings 300000; 
done