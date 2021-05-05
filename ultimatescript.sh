#!/bin/bash -x

#$readarray -t data < "$DATA"
#$echo ${data[*]}

#BARTIMAEUS
ROUND="4"
START=/home/irene/Documents/Studium/PhD/Poster_and_Praesentation/2020_Poster_EMBL/ernwin
NAME="CopStems"
PROGS=/home/irene/Programs/GitHub/RNA-Interaction-Workflow
SIMRNA=/home/irene/Programs/SimRNA_64bitIntel_Linux
DATA=/home/irene/Programs/SimRNA_64bitIntel_Linux/data
WHERE="local"
BUFFER="2"

#CORIDAN
#ROUND="3"
#START=/scratch/irene/Data_interaction/TESTFINAL
#NAME="1zci"
#PROGS=/scratch/irene/Data_interaction/RNA-Interaction-Workflow
#SIMRNA=/home/mescalin/irene/Programs/SimRNA_64bitIntel_Linux
#DATA=/home/mescalin/irene/Programs/SimRNA_64bitIntel_Linux/data
#WHERE="local"
#BUFFER="2"

SEQ="${START}/${NAME}.seq"
SS0="${START}/${NAME}.ss"
basename "$SS0"
BASESS0="$(basename -- $SS0)"

echo $ROUND
echo $START
echo $NAME
echo $PROGS
echo $SIMRNA
echo $DATA
echo $SEQ
echo $SS0

cp -r $DATA .


#####EXPAND#####
CURRENTROUND='3'

echo $CURRENTROUND
ROLD="$(($CURRENTROUND-1))"
NAMESS=$START/"${NAME}_${CURRENTROUND}.ss"
python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d "${START}/${ROLD}/${NAME}_${ROLD}.ss" -o $NAMESS -r -l

if [ $CURRENTROUND -eq 1 ]; then
	diff -q $NAMESS $SS0 1>/dev/null
else
	diff -q $NAMESS "${START}_${ROLD}.ss" 1>/dev/null
fi
if [[ $? == "0" ]]; then
		echo "The same"
	rm $NAMESS
	break 1
else
	echo "Not the same = continue"
	mkdir ${CURRENTROUND}
	mv $NAMESS ${START}/${CURRENTROUND}/
fi

$PROGS/SimRNA_scripts/job_simrna_start_expand.sh $START/${CURRENTROUND}/ $NAME $CURRENTROUND $SIMRNA expand $WHERE step
wait

mkdir $START/${CURRENTROUND}/analyse
for step in {1..10}; do
	mkdir ${START}/${CURRENTROUND}/analyse/${step}
	TRAFL="${NAME}_expand_${CURRENTROUND}_${step}_${step}.trafl"
	PDB="${NAME}_expand_${CURRENTROUND}_${step}_${step}-000001.pdb"
	cp $START/${CURRENTROUND}/${TRAFL} ${START}/${CURRENTROUND}/analyse/${step}/
	cp ${START}/${CURRENTROUND}/${PDB} ${START}/${CURRENTROUND}/analyse/${step}/
	SimRNA_trafl2pdbs $START/${CURRENTROUND}/analyse/${step}/*.pdb $START/${CURRENTROUND}/analyse/${step}/*.trafl :

	if [ $step = "1" ]; then
		python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_expand_${CURRENTROUND}_${step}.csv -u ${NAME}_expand_${CURRENTROUND}.csv -m 'w' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
	else
		python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_expand_${CURRENTROUND}_${step}.csv -u ${NAME}_expand_${CURRENTROUND}.csv -m 'a' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
	fi
done

mv ${NAME}_expand* ${START}/${CURRENTROUND}/

python ${PROGS}/continoussearch.py -p ${START}/${CURRENTROUND}/ --print --first "${NAME}_${CURRENTROUND}.ss" --second "${NAME}_${CURRENTROUND}.ss_cc" --interaction -i "${NAME}_expand_${CURRENTROUND}"

SSCCBP=(${START}/${CURRENTROUND}/*.ss_detected.bp)
BASIS=(${SSCCBP##*/})
NR="$(cut -d'-' -f2 <<<${BASIS})"
NR="$(cut -d'.' -f1 <<<${NR})"
ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
echo $BASIS
echo $NR
echo $ROUNDSEL

mv ${NAME}_expand* ${START}/${CURRENTROUND}/

SimRNA_trafl2pdbs $START/${CURRENTROUND}/"${NAME}_expand_${CURRENTROUND}_${CURRENTROUND}-000001.pdb" $START/${CURRENTROUND}/"${NAME}_expand_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA

ROUNDNEW="$(($CURRENTROUND+"1"))"
mkdir $START/$ROUNDNEW

SSCC=${START}/${CURRENTROUND}/"${NAME}_expand_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected"
cp $SSCC ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.ss_cc"

NEWPDB=${START}/${CURRENTROUND}/"${NAME}_expand_${CURRENTROUND}_**_AA.pdb"
cp $NEWPDB ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.pdb"

rm -r "${START}/${CURRENTROUND}/analyse"


#rm -r data
