#!/bin/bash -x

#$readarray -t data < "$DATA"
#$echo ${data[*]}

#BARTIMAEUS
#ROUND="2"
#START=/home/irene/Documents/TEST
#NAME="1zci"
#PROGS=/home/irene/Programs/GitHub/RNA-Interaction-Workflow
#SIMRNA=/home/irene/Programs/SimRNA_64bitIntel_Linux
#DATA=/home/irene/Programs/SimRNA_64bitIntel_Linux/data
#WHERE="cluster"


#CORIDAN
ROUND="3"
START=/scratch/irene/Data_interaction/TESTFINAL
NAME="1zci"
PROGS=/scratch/irene/Data_interaction/RNA-Interaction-Workflow
SIMRNA=/home/mescalin/irene/Programs/SimRNA_64bitIntel_Linux
DATA=/home/mescalin/irene/Programs/SimRNA_64bitIntel_Linux/data
WHERE="local"
BUFFER="2"

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

#"${PROGS}/SimRNA_scripts/job_simrna_start_initial.sh" $START $NAME "0" $SIMRNA "initial" $WHERE "step"
wait

cat $START/*.trafl > $START/"${NAME}_0_cat.trafl"
python "${PROGS}/traflminE.py" -i $START/"${NAME}_0_cat.trafl" -p $START/ -o "${NAME}_0_minE.trafl"
#trafl_extract_lowestE_frame.py "${START}/${NAME}_0_cat.trafl" #python2
rm $START/"${NAME}_0_cat.trafl"
SimRNA_trafl2pdbs $START/"${NAME}_initial_0_1_1-000001.pdb" $START/"${NAME}_0_minE.trafl" : AA

SSCC=${START}/"${NAME}_0_minE-000001.ss_detected"
cp $SSCC ${START}/"${NAME}_0.ss_cc"
SSCC=${START}/"${NAME}_0.ss_cc"
basename "$SSCC"
BASESSCC="$(basename -- $SSCC)"

NEWPDB=${START}/"${NAME}_0_minE-000001_AA.pdb"
cp $NEWPDB ${START}/"${NAME}_0.pdb"

####SURFACE FROM INITIAL####
$PROGS/SimRNA_scripts/job_simrna_start_surface.sh $START $NAME 0 $SIMRNA surface $WHERE step
wait

mkdir $START/analyse
for step in {1..10}; do
	mkdir $START/analyse/$step
	TRAFL="${NAME}_surface_0_${step}_${step}.trafl"
	PDB="${NAME}_surface_0_${step}_${step}-000001.pdb"
	cp $START/$TRAFL $START/analyse/${step}/
	cp $START/$PDB $START/analyse/${step}/
	SimRNA_trafl2pdbs $START/analyse/${step}/*.pdb $START/analyse/${step}/*.trafl :

	if [ $step = "1" ]; then
		python $PROGS/SSalignment.py -p $START/analyse/${step}/ -i $SSCC -c $SS0 -o "${NAME}_surface_0_${step}.csv" -u "${NAME}_surface_0.csv" -m 'w' -t $START/analyse/${step}/${TRAFL}
	else
		python $PROGS/SSalignment.py -p $START/analyse/${step}/ -i $SSCC -c $SS0 -o "${NAME}_surface_0_${step}.csv" -u "${NAME}_surface_0.csv" -m 'a' -t $START/analyse/${step}/${TRAFL}
	fi
done
mv ${NAME}_surface_0* ${START}/

python ${PROGS}/continoussearch.py -p ${START} --print --first $BASESS0 --second $BASESSCC --interaction -i "${NAME}_surface_0"

SSCCBP=(${START}/*.ss_detected.bp)
BASIS=(${SSCCBP##*/})
NR="$(cut -d'-' -f2 <<<${BASIS})"
NR="$(cut -d'.' -f1 <<<${NR})"
ROUNDSEL="$(cut -d'_' -f4 <<<${BASIS})"
echo $BASIS
echo $NR
echo $ROUNDSEL

SimRNA_trafl2pdbs $START/${NAME}_surface_0_${ROUNDSEL}_${ROUNDSEL}-000001.pdb $START/${NAME}_surface_0_${ROUNDSEL}_${ROUNDSEL}.trafl ${NR} AA

mv ${NAME}_surface_0* ${START}/

mkdir $START/1
SSCC=${START}/"${NAME}_surface_0_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected"
cp $SSCC ${START}/1/"${NAME}_1.ss_cc"
NEWPDB=${START}/"${NAME}_surface_0_**_AA.pdb"
cp $NEWPDB ${START}/1/"${NAME}_1.pdb"

rm -r "${START}/analyse"

#####EXPAND#####
for CURRENTROUND in `seq 1 1 ${ROUND}`; do
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
		mv $NAMESS ${START}/${CURRENTROUND}/
	fi

	$PROGS/SimRNA_scripts/job_simrna_start_expand.sh $START/${CURRENTROUND}/ $NAME $CURRENTROUND $SIMRNA expand_long00 $WHERE step
	wait

	mkdir $START/${CURRENTROUND}/analyse
	for step in {1..10}; do
		mkdir ${START}/${CURRENTROUND}/analyse/${step}
		TRAFL="${NAME}_expand_long00_${CURRENTROUND}_${step}_${step}.trafl"
		PDB="${NAME}_expand_long00_${CURRENTROUND}_${step}_${step}-000001.pdb"
		cp $START/${CURRENTROUND}/${TRAFL} ${START}/${CURRENTROUND}/analyse/${step}/
		cp ${START}/${CURRENTROUND}/${PDB} ${START}/${CURRENTROUND}/analyse/${step}/
		SimRNA_trafl2pdbs $START/${CURRENTROUND}/analyse/${step}/*.pdb $START/${CURRENTROUND}/analyse/${step}/*.trafl :

		if [ $step = "1" ]; then
			python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_expand_long00_${CURRENTROUND}_${step}.csv -u ${NAME}_expand_long00_${CURRENTROUND}.csv -m 'w' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
		else
			python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_expand_long00_${CURRENTROUND}_${step}.csv -u ${NAME}_expand_long00_${CURRENTROUND}.csv -m 'a' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
		fi
	done

	mv ${NAME}_expand_long00* ${START}/${CURRENTROUND}/
	
	python ${PROGS}/continoussearch.py -p ${START}/${CURRENTROUND}/ --print --first "${NAME}_${CURRENTROUND}.ss" --second "${NAME}_${CURRENTROUND}.ss_cc" --interaction -i "${NAME}_expand_long00_${CURRENTROUND}"

	SSCCBP=(${START}/${CURRENTROUND}/*.ss_detected.bp)
	BASIS=(${SSCCBP##*/})
	NR="$(cut -d'-' -f2 <<<${BASIS})"
	NR="$(cut -d'.' -f1 <<<${NR})"
	ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
	echo $BASIS
	echo $NR
	echo $ROUNDSEL

	mv ${NAME}_expand_long00* ${START}/${CURRENTROUND}/

	SimRNA_trafl2pdbs $START/${CURRENTROUND}/"${NAME}_expand_long00_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/${CURRENTROUND}/"${NAME}_expand_long00_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA

	ROUNDNEW="$(($CURRENTROUND+"1"))"
	mkdir $START/$ROUNDNEW

	SSCC=${START}/${CURRENTROUND}/"${NAME}_expand_long00_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected"
	cp $SSCC ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.ss_cc"

	NEWPDB=${START}/${CURRENTROUND}/"${NAME}_expand_long00_${CURRENTROUND}_**_AA.pdb"
	cp $NEWPDB ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.pdb"

	rm -r "${START}/${CURRENTROUND}/analyse"
done

rm -r data
