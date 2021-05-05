#!/bin/bash

START=/home/irene/Programs/GitLab/bvdv-pk/BVDVsegment_180-204_357-425/SimRNA/RE3/clustering/angle-measurement/set1
#START=/home/irene/Programs/GitLab/bvdv-pk/BVDVsegment_180-204_357-425/SimRNA/RE3/clustering/angle-measurement/set2
#START=/home/irene/Programs/GitLab/bvdv-pk/BVDVsegment_180-204_357-425/SimRNA/RE3/clustering/angle-measurement/set3
#START=/home/irene/Programs/GitLab/bvdv-pk/BVDVsegment_180-204_357-425/SimRNA/RE3/clustering/angle-measurement/set4
START=/home/irene/Programs/GitLab/bvdv-pk/BVDVsegment_180-204_357-425/SimRNA/RE6/clustering/todd_thrs5.00A_clust01
FORGI=/home/irene/Programs/GitHub/forgi/examples/
#NAME="todd_thrs5.00A"

echo $START
echo $NAME

#for pdb in $START/*.pdb; do
for pdb in $START/*_AA.pdb; do
	echo $pdb

	BASIS=(${pdb##*/})
	BASIS="${BASIS%.*}"
	echo $BASIS

	rnaConvert.py $pdb --pseudoknots -T fasta --to-file
	rnaConvert.py $pdb --pseudoknots -T forgi --to-file
	rnaConvert.py $pdb --pseudoknots -T element_string --to-file
	mv ${BASIS}* ${START}/

	x3dna-dssr --input=${pdb} --json -o="${START}/${BASIS}.json"

	python ${FORGI}/coaxial_stacking.py "${START}/${BASIS}_A.cg" --dssr-json "${START}/${BASIS}.json" -l --csv "${START}/set1.cvs"
	python ${FORGI}/interior_loop_angles.py ${pdb} -o "${START}/${BASIS}_intloop.cvs"

done
