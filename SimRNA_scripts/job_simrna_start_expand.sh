#!/bin/bash

#iterations 10,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand local|cluster |random

#iterations 5,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long00 local|cluster |random

#iterations 20,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long01 local|cluster |random

#iterations 30,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long02 local|cluster |random

#iterations 100,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long03 local|cluster |random

#iterations 300,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long04 local|cluster |random

#iterations 600,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long05 local|cluster |random

#iterations 1,000,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long06 local|cluster |random

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5
WHERE=$6

PDB="${INOUTDIR}/${INITNAME}_${ROUND}.pdb"
echo "$PDB" "$INITNAME" "$ROUND"
SSFILE="${INOUTDIR}/${INITNAME}_${ROUND}.ss"
SIMRNACONFIG=${SIMRNA}/"config_"${CONFIG}".dat"
SIMRNADATA="${SIMRNA}/data"
NAME="${INITNAME}_${CONFIG}_${ROUND}"

# checkout the right config file
IFS='_' read -r -a pos <<< "$CONFIG"
#START="/scratch/irene/Data_interaction/RNA-Interaction-Workflow/SimRNA_scripts/job_simrna_script_"$pos".sh"
#Bartimaeus:
START="/home/irene/Programs/GitHub/RNA-Interaction-Workflow/SimRNA_scripts/job_simrna_script_"$pos".sh"

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
if [ "$WHERE" = "local" ]; then
	ln -s "$SIMRNADATA" .
	if [ "$7" = "random" ]; then
		for step in {1..10}; do
			random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
			# Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${random}"
			SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$random" -o "$NEWNAME" >& "$NEWNAME".log &
		done
	else
		for step in {1..10}; do
		#for step in {1..2}; do #Scenario testing
			NEWNAME="${INITNAME}_${CONFIG}_${ROUND}_${step}_${step}"
			#echo SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$step" -o "$NEWNAME" >& "$NEWNAME".log &
			SimRNA -p "$PDB" -S "$SSFILE" -c "$SIMRNACONFIG" -R "$step" -o "$NEWNAME" >& "$NEWNAME".log &
		done
	fi
	wait
	#rm -r data
	ls "$NAME"*
	mv "$NAME"* "$INOUTDIR"

elif [ "$WHERE" = "cluster" ]; then
	sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$START" "$PDB" "$SSFILE" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" "$7"
	wait
else
	echo "No calculation location given local|cluster"
fi
wait
