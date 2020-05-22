#!/bin/bash

#iterations 10,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand

#iterations 20,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long01
 
#iterations 30,000 stepsize 1
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long02
 
#iterations 100,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long03 

#iterations 300,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long04 

#iterations 600,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long05 

#iterations 1,000,000 stepsize 100
##SimRNA_scripts/job_simrna_start_expand.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ expand_long06 

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5

NAME="${INITNAME}_${CONFIG}_${ROUND}"
PDB="${INOUTDIR}/${INITNAME}.pdb"
SSFILE="${INOUTDIR}/${INITNAME}_${ROUND}.ss"
SIMRNACONFIG=${SIMRNA}/"config_"${CONFIG}.dat
SIMRNADATA=${SIMRNA}/data/

IFS='_' read -r -a pos <<< "$CONFIG"
START="SimRNA_scripts/job_simrna_script_"$pos".sh"

#Check if directories/files exist and all filenames/directories are correct
echo 'Outputname' "$NAME"

if [ ! -f "$START" ]; then
    echo 'Sequence file does not exist!' "$START"
    exit
    else
    echo "$START"
fi

if [ ! -d "$INOUTDIR" ]; then
    echo 'In-/Output directory does not exist!' "$INOUTDIR"
    exit
    else
    echo "$INOUTDIR"
fi

if [ ! -d "$SIMRNADATA" ]; then
    echo 'SimRNA data directory does not exist!' "$SIMRNADATA" 
    exit
    else
    echo "$SIMRNADATA"
fi

if [ ! -f "$PDB" ]; then
    echo 'PDB file does not exist!' "$PDB"
    exit
    else
    echo "$PDB"
fi

if [ ! -f "$SSFILE" ]; then
    echo 'Secondary structure file does not exist!' "$SSFILE"
    exit
    else
    echo "$SSFILE"
fi

if [ ! -f "$SIMRNACONFIG" ]; then
    echo 'SimRNA config file does not exist!' "$SIMRNACONFIG"
    exit
    else
    echo "$SIMRNACONFIG"
fi

  sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$START" "$PDB" "$SSFILE" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" &

wait

