#!/bin/bash

##SimRNA_scripts/job_simrna_start_surface.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ surface 

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5

NAME="${INITNAME}_${CONFIG}_${ROUND}"
PDB="${INOUTDIR}/${INITNAME}_${ROUND}.pdb"
SIMRNACONFIG=${SIMRNA}/"config_"${CONFIG}.dat
SIMRNADATA=${SIMRNA}/data/

START="SimRNA_scripts/job_simrna_script_"${CONFIG}".sh"

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

if [ ! -f "$SIMRNACONFIG" ]; then
    echo 'SimRNA config file does not exist!' "$SIMRNACONFIG"
    exit
    else
    echo "$SIMRNACONFIG"
fi

  sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$START" "$PDB" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" &

wait