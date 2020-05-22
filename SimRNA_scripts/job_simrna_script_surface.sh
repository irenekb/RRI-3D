#!/bin/bash
#SBATCH --job-name=SimRNAsurface
#SBATCH --time=144:0:0
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-10

#Get commandline arguments:
#CSV1: PDB
#CSV2: SimRNA Data conection, CSV3: SimRNA configfile
#CSV4: OUTDIR, CSV5: NAME
PDB="$1"
CSV2="$2"
CSV3="$3"
OUTPUT_DIRECTORY="$4"
NAME="$5"

NAME+="_"${SLURM_ARRAY_TASK_ID}

#Change to local directory
cd "$TMPDIR"

#Create link/copy file to file (SimRNA needs it!!!)
cp -r "$CSV2" .

#Run SimRNA
echo SimRNA -p "$PDB" -c "$CSV3" -R ${SLURM_ARRAY_TASK_ID} -o "$NAME" >& "$NAME".log 

#Copy files to output directory
rm -r data
cp * "$OUTPUT_DIRECTORY"