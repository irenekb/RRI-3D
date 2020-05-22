#!/bin/bash
#SBATCH --job-name=SimRNAexpand
#BATCH --time=144:0:0
#BATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-10

#Get commandline arguments:
#CSV1:pdb, CSV2: secondary-structure file 
#CSV3: SimRNA Data conection, CSV4: SimRNA configfile
#CSV5: OUTDIR, CSV6: NAME
PDB="$1"
SSFILE="$2"
CSV3="$3"
CSV4="$4"
OUTPUT_DIRECTORY="$5"
NAME="$6"

NAME+="_"${SLURM_ARRAY_TASK_ID}

#Change to local directory
cd "$TMPDIR"

#Create link/copy file to file (SimRNA needs it!!!)
cp -r "$CSV3" .

#Run SimRNA
echo SimRNA -p "$PDB" -S "$SSFILE" -c "$CSV4" -R ${SLURM_ARRAY_TASK_ID} -o "$NAME" >& "$NAME".log 

#Copy files to output directory
rm -r data
cp * "$OUTPUT_DIRECTORY"
