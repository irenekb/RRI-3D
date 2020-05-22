#!/bin/bash

##SimRNA_scripts/job_simrna_start_initial.sh SimRNA_interaction/1zci/ 1zci 00 ~/Programs/SimRNA_64bitIntel_Linux/ initial 

INOUTDIR=$(realpath "$1")
INITNAME="$2"
ROUND="$3"
SIMRNA=$(realpath "$4")
CONFIG=$5

NAME="${INITNAME}_${CONFIG}_${ROUND}"
SEQFILE="${INOUTDIR}/${INITNAME}.seq"
SSFILE="${INOUTDIR}/${INITNAME}_${ROUND}.ss"
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

if [ ! -f "$SEQFILE" ]; then
    echo 'Sequence file does not exist!' "$SEQFILE"
    exit
    else
    echo "$SEQFILE"
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

  sbatch --output "$INOUTDIR"/"$NAME"_%j.log "$START" "$SEQFILE" "$SSFILE" "$SIMRNADATA" "$SIMRNACONFIG" "$INOUTDIR" "$NAME" &

wait

