#!/bin/bash

#SBATCH --job-name=SimRNAsurface
#SBATCH --time=00:20:00
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=6

## script started from job_simrna_start_surface.sh for a cluster calculation
#Get commandline arguments:
#CSV1: PDB
#CSV2: SimRNA Data conection, CSV3: SimRNA configfile
#CSV4: OUTDIR, CSV5: NAME
PDB="$1"
CSV2="$2"
CSVCONFIG3="$3"
OUTPUT_DIRECTORY="$4"
SEEDTYPE="$6"

#Change to local directory
cd "$TMPDIR"

#Create link/copy file to file (SimRNA needs it!!!)
cp -r "$CSV2" .

#Run SimRNA
if [ "$SEEDTYPE" = "random" ]; then
		random=$(od -N 4 -t uL -An < /dev/urandom | tr -d " ")
    # Reads 4 bytes from the random device and formats them as unsigned integer between 0 and 2^32-1
		NAME="$5_${SLURM_ARRAY_TASK_ID}_${random}"
		SimRNA -p "$PDB"  -c "$CONFIG" -R "$random" -o "$NAME" >& "$NAME".log
	done
else
		NAME="$5_${SLURM_ARRAY_TASK_ID}_${SLURM_ARRAY_TASK_ID}"
		SimRNA -p "$PDB" -c "$CONFIG" -R ${SLURM_ARRAY_TASK_ID} -o "$NAME" >& "$NAME".log
	done
fi


#Copy files to output directory
rm -r data
cp * "$OUTPUT_DIRECTORY"
