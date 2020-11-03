#!/bin/bash -x

#$readarray -t data < "$DATA"
#$echo ${data[*]}

#BARTIMAEUS

START=/scr/coridan/irene/Data_interaction/SimRNA/CopStems/CopStemsErnwin/best
NAME="CopStemsErnwin900"
EXTEND="expand"
PROGS=/scr/coridan/irene/Data_interaction/RNA-Interaction-Workflow
SIMRNA=~/Programs/SimRNA_64bitIntel_Linux/
DATA=~/Programs/SimRNA_64bitIntel_Linux/data

SEQ="${START}/${NAME}.seq"
SS0="${START}/${NAME}_00.ss"
basename "$SS0"
BASESS0="$(basename -- $SS0)"

echo $START
echo $NAME
echo $PROGS
echo $SIMRNA
echo $DATA
echo $SEQ
echo $SS0

ln -s $DATA .


mkdir $START/${NAME}

NAMESS=$START/"${NAME}_00.ss"
NAMECC=$START/"${NAME}_00.ss_cc"

#$PROGS/SimRNA_scripts/job_simrna_start_expand.sh $START/00/ $NAME $CURRENTROUND $SIMRNA ${EXTEND} $WHERE step
#wait
mkdir $START/$NAME/
mkdir $START/$NAME/analyse
for step in {1..10}; do
	mkdir ${START}/$NAME/analyse/${step}
	TRAFL="${NAME}_${EXTEND}_00_${step}_${step}.trafl"
	PDB="${NAME}_${EXTEND}_00_${step}_${step}-000001.pdb"
	cp $START/${TRAFL} ${START}/$NAME/analyse/${step}/
	cp ${START}/${PDB} ${START}/$NAME/analyse/${step}/
	SimRNA_trafl2pdbs $START/$NAME/analyse/${step}/*.pdb $START/$NAME/analyse/${step}/*.trafl :

	if [ $step = "1" ]; then
		python $PROGS/SSalignment.py -p $START/$NAME/analyse/${step}/ -i $NAMESS -c $NAMESS -o ${NAME}_${EXTEND}_00_${step}.csv -u ${NAME}_${EXTEND}_00.csv -m 'w' -t $START/$NAME/analyse/${step}/$TRAFL
	else
		python $PROGS/SSalignment.py -p $START/$NAME/analyse/${step}/ -i $NAMESS -c $NAMESS -o ${NAME}_${EXTEND}_00_${step}.csv -u ${NAME}_${EXTEND}_00.csv -m 'a' -t $START/$NAME/analyse/${step}/$TRAFL
	fi
done

mv ${NAME}_${EXTEND}* ${START}/${NAME}/

python ${PROGS}/continoussearch.py -p ${START}/${NAME}/ --print --first $BASESS0 --interaction -i "${NAME}_${EXTEND}_00"

SSCCBP=(${START}/${NAME}/*.ss_detected.bp)
BASIS=(${SSCCBP##*/})
NR="$(cut -d'-' -f2 <<<${BASIS})"
NR="$(cut -d'.' -f1 <<<${NR})"
ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
echo $BASIS
echo $NR
echo $ROUNDSEL

mv ${NAME}_${EXTEND}* ${START}/$NAME

SimRNA_trafl2pdbs $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA

mv $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}-${NR}.pdb" $START/$NAME/

rm -r data





