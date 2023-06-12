#!/bin/bash

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5
WHERE=$6

NAME="${INITNAME}_${CONFIG}_${ROUND}"
SEQFILE="${INOUTDIR}/${INITNAME}.seq"
SSFILE="${INOUTDIR}/${INITNAME}_${ROUND}.ss"
SIMRNACONFIG=${SIMRNA}/"config_"${CONFIG}".dat"
SIMRNADATA=${SIMRNA}/data/

START="/scratch/irene/Data_interaction/RNA-Interaction-Workflow/SimRNA_scripts/job_simrna_script_"${CONFIG}".sh"

#Check if directories/files exist and all filenames/directories are correct
errors=()

if [ ! -f "$START" ]; then
	errors+=('Starting script does not exist!' "$START")
else
       	echo "$START"
fi

if [ ! -d "$INOUTDIR" ]; then
	errors+=('In-/Output directory does not exist!' "$INOUTDIR")
else
       	echo "$INOUTDIR"
fi

if [ ! -d "$SIMRNADATA" ]; then
	errors+=('SimRNA data directory does not exist!' "$SIMRNADATA") 
else
	echo "$SIMRNADATA"
fi

if [ ! -f "$SEQFILE" ]; then
	errors+=('Sequence file does not exist!' "$SEQFILE")
else
	echo "$SEQFILE"
fi

if [ ! -f "$SSFILE" ]; then
	errors+=('Secondary structure file does not exist!' "$SSFILE")
else
	echo "$SSFILE"
fi

if [ ! -f "$SIMRNACONFIG" ]; then
	errors+=('SimRNA config file does not exist!' "$SIMRNACONFIG")
else
	echo "$SIMRNACONFIG"
fi

 
# if there is an error print it and exit
if [ ${#errors[@]} -eq 0 ] ; then
	echo 'Outputname' "$NAME"
else
	shuf -e "${errors[@]}"
	exit
fi


# start the job  localy or on the cluster:
if [ "$WHERE" = "local" ]; then
	#ln -s "$SIMRNADATA" .
	if [ "$7" = "random" ]; then
		for step in {1..10}; do
			random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
			# Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1	
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${random}"
			SimRNA -s "$SEQFILE" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$random" -o "$NEWNAME" >& "$NEWNAME".log & 
		done
	else
		for step in {1..10}; do
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${step}"
			SimRNA -s "$SEQFILE" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$step" -o "$NEWNAME" >& "$NEWNAME".log & 
		done
	fi
	wait
	#rm -r data
	ls "$NAME"*
	mv "$NAME"* "$INOUTDIR"

elif [ "$WHERE" = "cluster" ]; then
	sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$START" "$SEQFILE" "$SSFILE" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" "$7"
	wait
else
	echo "No calculation location given local|cluster"
fi
wait

