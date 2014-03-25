if [ -z $2 ]
then
	echo "usage:$0 <path to assembled files> <master assembly file> <crumbs>"
	exit
fi


ASSEMBLED_PATH=$1
OUTFILE=$2
CRUMBS=$3

# clear out any existing outfiles.
if [ -f $OUTFILE ]; then rm $OUTFILE; fi
    
assembled_file_list=`ls -1 $ASSEMBLED_PATH/*.assembled`
punted_file_list=`ls -1 $ASSEMBLED_PATH/*.punted`

for assembled_file in $assembled_file_list
do
    echo "adding" $assembled_file "to master SV assembly file ("$OUTFILE")"
    awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,\
                          $11,$12,$13,$14,$15,$16,$17,$18,\
                          $19,FILENAME"_"$20,$21,$22,$23}' $assembled_file \
    >> $OUTFILE

done

if [[ -z "$3" ]]; then
	echo -e  "\tCleaning up old files..."
	for punted_file in $punted_file_list
	do
    		`rm $punted_file`

    	done
fi

if [[ -z "$3" ]]; then
	echo -e  "\tCleaning up old files..."
	for assembled_file in $assembled_file_list
	do
		`rm $assembled_file`

    	done
fi

