#!/bin/bash

#SBATCH --job-name=SimRNAexpand
#BATCH --time=01:0:0
#BATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10

## script started from job_simrna_start_expand.sh for a cluster calculation
#Get commandline arguments:
#CSV1:pdb, CSV2: secondary-structure file 
#CSV3: SimRNA Data conection, CSV4: SimRNA configfile
#CSV5: OUTDIR, CSV6: NAME
PDB="$1"
SSFILE="$2"
DATA="$3"
CONFIG="$4"
OUTPUT_DIRECTORY="$5"
SEEDTYPE="$7"


#Change to local directory
cd "$TMPDIR"


#Create link/copy file to file (SimRNA needs it!!!)
cp -r "$DATA" .


#Run SimRNA
if [ "$SEEDTYPE" = "random" ]; then
        for step in {1..10}; do
		random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
                # Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1
		NAME="$6_${step}_${random}"
		SimRNA -p "$PDB" -S "$SSFILE" -c "$CONFIG" -R "$random" -o "$NAME" >& "$NAME".log & 
	done
else
	for step in {1..10}; do	
		NAME="$6_${step}_${step}"
		SimRNA -p "$PDB" -S "$SSFILE" -c "$CONFIG" -R "$step" -o "$NAME" >& "$NAME".log &
	done
fi


#Copy files to output directory
rm -r data
cp * "$OUTPUT_DIRECTORY"
