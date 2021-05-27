#!/bin/bash -x

#$readarray -t data < "$DATA"
#$echo ${data[*]}

#BARTIMAEUS
#ROUND="0"
##START=/home/irene/Documents/Studium/PhD/Poster_and_Praesentation/2020_Poster_EMBL/ernwin
#START=/home/irene/Programs/GitHub/RNA-Interaction-Workflow/CALCULATION_NO-UPDATE/CopStems/simrna
#NAME="CopStems"
#PROGS=/home/irene/Programs/GitHub/RNA-Interaction-Workflow
#SIMRNA=/home/irene/Programs/SimRNA_64bitIntel_Linux
#DATA=/home/irene/Programs/SimRNA_64bitIntel_Linux/data
#WHERE="local"
#BUFFER="2"
#TOTAL="NO"


TEST
ROUND="0"
ROUNDS="100"
#START=/home/irene/Programs/GitHub/RNA-Interaction-Workflow/CALCULATION_NO-UPDATE/CopStems/simrna
START=/home/irene/Programs/GitHub/RNA-Interaction-Workflow/CALCULATION_NO-UPDATE/1CZI/simrna
NAME="1zci"
PROGS=/home/irene/Programs/GitHub/RNA-Interaction-Workflow
SIMRNA=/home/irene/Programs/SimRNA_64bitIntel_Linux
DATA=/home/irene/Programs/SimRNA_64bitIntel_Linux/data
WHERE="local"
BUFFER="2"
TOTAL="NO"
TYPE="expand"
EXTEND="expand_long00" #expand_long01

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
PDB="${START}/${NAME}_${ROUND}.pdb"
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

ln -s $DATA .

NAMESS=$START/"${NAME}_${ROUND}.ss"
NAMECC=$START/"${NAME}_${ROUND}.ss_cc"



#####START from ernwin reconstructed pdb#####
echo "$PROGS/SimRNA_scripts/job_simrna_start_${TYPE}.sh" $START $NAME $ROUND $SIMRNA $EXTEND $WHERE "step"
"$PROGS/SimRNA_scripts/job_simrna_start_${TYPE}.sh" $START $NAME $ROUND $SIMRNA $EXTEND $WHERE "step"
wait

mkdir $START/$ROUND
mkdir $START/$ROUND/analyse
for step in {1..10}; do
#for step in {1..2}; do #Scenario testing
	mkdir ${START}/$ROUND/analyse/${step}

	TRAFL="${NAME}_${EXTEND}_${ROUND}_${step}_${step}.trafl"
	PDB="${NAME}_${EXTEND}_${ROUND}_${step}_${step}-000001.pdb"
	cp $START/${TRAFL} ${START}/$ROUND/analyse/${step}/
	cp ${START}/${PDB} ${START}/$ROUND/analyse/${step}/
	SimRNA_trafl2pdbs $START/$ROUND/analyse/${step}/*.pdb $START/$ROUND/analyse/${step}/*.trafl :

	if [ $step = "1" ]; then
		python $PROGS/SSalignment.py -p $START/$ROUND/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${EXTEND}_0_${step}.csv -u ${NAME}_${EXTEND}_0.csv -m 'w' -t $START/$ROUND/analyse/${step}/$TRAFL
	else
		python $PROGS/SSalignment.py -p $START/$ROUND/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${EXTEND}_0_${step}.csv -u ${NAME}_${EXTEND}_0.csv -m 'a' -t $START/$ROUND/analyse/${step}/$TRAFL
	fi
done

mv ${NAME}_${EXTEND}* ${START}/$ROUND

python ${PROGS}/continoussearch.py -p ${START}/$ROUND --print --first "${NAME}_${ROUND}.ss" --second "${NAME}_${ROUND}.ss_cc" -i "${NAME}_${EXTEND}_0" --force
#get the forced structure and search within the individual runs

SSCCBP=(${START}/${NAME}/*.ss_detected.bp)
BASIS=(${SSCCBP##*/})
NR="$(cut -d'-' -f2 <<<${BASIS})"
NR="$(cut -d'.' -f1 <<<${NR})"
ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
ROUNDSEL="$(cut -d'-' -f1 <<<${ROUNDSEL})"
echo $BASIS
echo $NR
echo $ROUNDSEL

SimRNA_trafl2pdbs $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA
ROUNDNEW="$(($ROUND+"1"))"
mkdir $START/${ROUNDNEW}

cp $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}_AA.pdb" $START/$ROUND/"${NAME}_${ROUNDNEW}.pdb"
cp $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected" $START/$ROUNDNEW/"${NAME}_${ROUNDNEW}.ss_cc"

#PDB result for the next round
ROLD="$(($ROUNDNEW-1))"
mv $START/$ROUND/"${NAME}_${ROUNDNEW}.pdb" $START/$ROUNDNEW/"${NAME}_${ROUNDNEW}.pdb"

# interaction lenght + bp
mv "${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.il" ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.il"
mv "${NAME}_${EXTEND}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.bp" ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.bp"

#clean up directory
mv "$START/${NAME}_${EXTEND}"* ${START}/${ROUND}/
mv "$START/${NAME}_0"* ${START}/${ROUND}/
#rm -r "${START}/${ROUND}/analyse"

#####EXPAND#####
for CURRENTROUND in `seq 1 1 ${ROUNDS}`; do
	ROLD="$(($CURRENTROUND-"1"))"
	echo $CURRENTROUND
	#expand SS
	NAMESS=$START/${CURRENTROUND}/"${NAME}_${CURRENTROUND}.ss"
	NAMESSOLD=$START/${ROLD}/"${NAME}_${ROLD}.ss"
	python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d "${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc" -o $NAMESS -r -l

	#Check if an expansion is possible
	if [ $CURRENTROUND -eq 1 ]; then
		diff -q $NAMESS $SS0 1>/dev/null
	else
		diff -q $NAMESS NAMESSOLD 1>/dev/null
	fi
	if [[ $? == "0" ]]; then
		echo "The same"
		rm -r data
		exit 0 		#break 1
	else
		echo "Not the same = continue"
	fi

	#start the SimRNA expansion job
	"$PROGS/SimRNA_scripts/job_simrna_start_${TYPE}.sh" $START/${CURRENTROUND}/ $NAME $CURRENTROUND $SIMRNA $EXTEND $WHERE step
	wait

	#analyse of the SimRNA run
	mkdir $START/${CURRENTROUND}/analyse
	for step in {1..10}; do
	#for step in {1..2}; do #Scenario testing
		mkdir ${START}/${CURRENTROUND}/analyse/${step}
		TRAFL="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}.trafl"
		PDB="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}-000001.pdb"
		cp $START/${CURRENTROUND}/${TRAFL} ${START}/${CURRENTROUND}/analyse/${step}/
		cp ${START}/${CURRENTROUND}/${PDB} ${START}/${CURRENTROUND}/analyse/${step}/
		SimRNA_trafl2pdbs $START/${CURRENTROUND}/analyse/${step}/*.pdb $START/${CURRENTROUND}/analyse/${step}/*.trafl :

		if [ $step = "1" ]; then
			#first round + create outputfiles
			python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'w' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
		else
			#append the following output
			python $PROGS/SSalignment.py -p $START/${CURRENTROUND}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'a' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL
		fi
	done

	mv ${NAME}_${EXTEND}* ${START}/${CURRENTROUND}/

	#find the starting structure for the next round
	python ${PROGS}/continoussearch.py -p ${START}/${CURRENTROUND} --print --first "${NAME}_${CURRENTROUND}.ss" --second "${NAME}_${CURRENTROUND}.ss_cc" --interaction -i "${NAME}_${EXTEND}_${CURRENTROUND}"


	SSCCBP=(${START}/${CURRENTROUND}/*.ss_detected.bp)
	BASIS=(${SSCCBP##*/})
	NR="$(cut -d'-' -f2 <<<${BASIS})"
	NR="$(cut -d'.' -f1 <<<${NR})"
	ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
	ROUNDSEL="$(cut -d'-' -f1 <<<${ROUNDSEL})"
	echo $BASIS
	echo $NR
	echo $ROUNDSEL

	#calculate a full atom strucutre
	SimRNA_trafl2pdbs $START/$CURRENTROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$CURRENTROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA


	#prepare everything for the next round
	ROUNDNEW="$(($CURRENTROUND+"1"))"
	mkdir $START/$ROUNDNEW

	NEWSSCC=${START}/${CURRENTROUND}/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected"
	cp $NEWSSCC ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.ss_cc"

	NEWPDB=$START/$CURRENTROUND/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}_AA.pdb"
	mv $NEWPDB $START/$ROUNDNEW/"${NAME}_${ROUNDNEW}.pdb"

	# interaction lenght + bp
	mv "${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.il" ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.il"
	mv "${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.bp" ${START}/${ROUNDNEW}/"${NAME}_${ROUNDNEW}.bp"

	#clean up directory
	#rm -r "${START}/${CURRENTROUND}/analyse"


	echo "Compare only interaction"
	declare -i ssjet=$(cat "${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.il" ) #typeset -i
	declare -i ssnew=$(cat "${START}/${ROUNDNEW}/${NAME}_${ROUNDNEW}.il")
	echo "Interaction length old: $ssold; Interaction length new: $ssjet"

	if [[ $ssnew -gt $ssjet ]]; then #greater than -gt
		echo "Not the same = continue"
	else
		echo "The same"
		rm -r data
		exit 0
	fi

done

echo "SIMULATION DONE"

rm -r data
