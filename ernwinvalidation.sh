#!/bin/bash -x

#$readarray -t data < "$DATA"
#$echo ${data[*]}

#BARTIMAEUS

#START=/scr/coridan/irene/Data_interaction/SimRNA/CopStems/CopStemsErnwin/best
START=/home/irene/Documents/Studium/PhD/m_virus/BVDV/translate-after-run/SimRNA_Start-BVDVsegment2pdb
NAME="BVDVsegment2"
EXTEND="surface"
ROUND="0"
WHERE="local"
#PROGS=/scr/coridan/irene/Data_interaction/RNA-Interaction-Workflow
#SIMRNA=~/Programs/SimRNA_64bitIntel_Linux/
#DATA=~/Programs/SimRNA_64bitIntel_Linux/data
PROGS=/home/irene/Programs/GitHub/RNA-Interaction-Workflow
SIMRNA=/home/irene/Programs/SimRNA_64bitIntel_Linux
DATA=/home/irene/Programs/SimRNA_64bitIntel_Linux/data


#SEQ="${START}/${NAME}.seq"
SEQ="${START}/${NAME}.pdb"
SS0="${START}/${NAME}.ss"
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


#mkdir $START/${NAME}

NAMESS=$START/"${NAME}_0.ss"
NAMECC=$START/"${NAME}_0.ss_cc"

#$PROGS/SimRNA_scripts/job_simrna_start_expand.sh $START/00/ $NAME $CURRENTROUND $SIMRNA ${EXTEND} $WHERE step
#wait
#echo $PROGS/SimRNA_scripts/job_simrna_start_surface.sh $START $NAME $ROUND $SIMRNA ${EXTEND} $WHERE
#$PROGS/SimRNA_scripts/job_simrna_start_surface.sh $START $NAME $ROUND $SIMRNA ${EXTEND} $WHERE
wait

#mkdir $START/$NAME/
#mkdir $START/$NAME/analyse
for step in {1..10}; do
	#mkdir ${START}/$NAME/analyse/${step}

	#TRAFL="${NAME}_${EXTEND}_00_${step}_${step}.trafl"
	#PDB="${NAME}_${EXTEND}_00_${step}_${step}-000001.pdb"

	TRAFL="${NAME}_${EXTEND}_0_${step}_${step}.trafl"
	PDB="${NAME}_${EXTEND}_0_${step}_${step}-000001.pdb"
	#cp $START/${TRAFL} ${START}/$NAME/analyse/${step}/
	#cp ${START}/${PDB} ${START}/$NAME/analyse/${step}/
	#SimRNA_trafl2pdbs $START/$NAME/analyse/${step}/*.pdb $START/$NAME/analyse/${step}/*.trafl :

	#if [ $step = "1" ]; then
		#python $PROGS/SSalignment.py -p $START/$NAME/analyse/${step}/ -i $NAMESS -c $NAMESS -o ${NAME}_${EXTEND}_0_${step}.csv -u ${NAME}_${EXTEND}_0.csv -m 'w' -t $START/$NAME/analyse/${step}/$TRAFL
	#else
		#python $PROGS/SSalignment.py -p $START/$NAME/analyse/${step}/ -i $NAMESS -c $NAMESS -o ${NAME}_${EXTEND}_0_${step}.csv -u ${NAME}_${EXTEND}_0.csv -m 'a' -t $START/$NAME/analyse/${step}/$TRAFL
	#fi
done

mv ${NAME}_${EXTEND}* ${START}/${NAME}/

#python ${PROGS}/continoussearch.py -p ${START}/${NAME}/ --print --first $BASESS0 --interaction -i "${NAME}_${EXTEND}_00"
python ${PROGS}/continoussearch.py -p ${START}/${NAME}/ --print --first $BASESS0 --interaction -i "${NAME}_${EXTEND}_0"
#python /home/irene/Programs/GitHub/RNA-Interaction-Workflow/continoussearch.py -p BVDVsegment2/ --print --first BVDVsegment2_0.ss --force -i BVDVsegment2_surface_0 --verbose

SSCCBP=(${START}/${NAME}/*.ss_detected.bp)
BASIS=(${SSCCBP##*/})
NR="$(cut -d'-' -f2 <<<${BASIS})"
NR="$(cut -d'.' -f1 <<<${NR})"
ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
echo $BASIS
echo $NR
echo $ROUNDSEL

mv ${NAME}_${EXTEND}* ${START}/$NAME

#SimRNA_trafl2pdbs $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA
#mv $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}-${NR}.pdb" $START/$NAME/

SimRNA_trafl2pdbs $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_00_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_0_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA
mv $START/$NAME/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_0_${ROUNDSEL}_${ROUNDSEL}-${NR}.pdb" $START/$NAME/


rm -r data
