#!/bin/bash

#SBATCH --job-name=SimRNAinitial
#SBATCH --time=24:0:0
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=10
#SBATCH --array=1-10

## script started from job_simrna_start_initial.sh for a cluster calculation

#Get commandline arguments:
#CSV1:sequence, CSV2: secondary-structure file,
#CSV3: SimRNA Data conection, CSV4: SimRNA configfile
#CSV5: OUTDIR, CSV6: NAME
SEQFILE="$1"
SSFILE="$2"
CSV3="$3"
CONFIG="$4"
OUTPUT_DIRECTORY="$5"
SEEDTYPE="$7"

#Change to local directory
cd "$TMPDIR"

#Create link/copy file to file (SimRNA needs it!!!)
cp -r "$CSV3" .

#Run SimRNA
if [ "$SEEDTYPE" = "random" ]; then
	random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
        # Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1
	NAME="$6_${SLURM_ARRAY_TASK_ID}_${random}"
	SimRNA -s "$SEQFILE" -S "$SSFILE" -c "$CONFIG" -R "$random" -o "$NAME" >& "$NAME".log
else
	NAME="$6_${SLURM_ARRAY_TASK_ID}_${SLURM_ARRAY_TASK_ID}"
	SimRNA -s "$SEQFILE" -S "$SSFILE" -c "$CONFIG" -R ${SLURM_ARRAY_TASK_ID} -o "$NAME" >& "$NAME".log
fi

#Copy files to output directory
rm -r data
cp * "$OUTPUT_DIRECTORY"
