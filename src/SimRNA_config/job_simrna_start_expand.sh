#!/bin/bash

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5
WHERE=$6
SEED=$7
SIMROUND=$8
SIMSTART=$9

PDB="${INOUTDIR}/${INITNAME}_${ROUND}.pdb"
SSFILE="${INOUTDIR}/${INITNAME}_${ROUND}.ss"
SIMRNACONFIG=${SIMRNA}/"config_"${CONFIG}".dat"
SIMRNADATA="${SIMRNA}/data"
NAME="${INITNAME}_${CONFIG}_${ROUND}"
echo "$PDB" "$INITNAME" "$ROUND" "$SSFILE" "$WHERE"

# checkout the right config file
IFS='_' read -r -a pos <<< "$CONFIG"
SIMSTART="${SIMSTART}/job_simrna_script_"$pos".sh"

#Check if directories/files exist and all filenames/directories are correct
errors=()

if [ ! -f "$SIMSTART" ]; then
	errors+=('Starting script does not exist!' "$SIMSTART")
else
	echo "$SIMSTART"
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

if [ ! -f "$PDB" ]; then
	errors+=('PDB file does not exist!' "$PDB")
else
	echo "$PDB"
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
if [ $WHERE = "local" ]; then
	ln -s "$SIMRNADATA" .
	if [ $SEED = "random" ]; then
		#for step in {1..10}; do
		for step in $(seq 1 ${SIMROUND}); do
			random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
			# Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${random}"
			echo "SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$random" -o "$NEWNAME" >& "$NEWNAME".log &"
			SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$random" -o "$NEWNAME" >& "$NEWNAME".log &
		done
	else
		for step in $(seq 1 ${SIMROUND}); do
		#for step in {1..10}; do
		#for step in {1..2}; do #Scenario testing
			echo $step
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${step}"
			echo "SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$step" -o "$NEWNAME" >& "$NEWNAME".log &"
			SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$step" -o "$NEWNAME" >& "$NEWNAME".log &
		done
	fi
	wait
	#rm -r data
	ls "$NAME"*
	mv "$NAME"* "$INOUTDIR"

elif [ $WHERE = "cluster" ]; then
	sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$SIMSTART" "$PDB" "$SSFILE" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" "$SEED" "$SIMROUND"
	wait
else
	echo "No calculation location given local|cluster"
fi
wait
