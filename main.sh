#!/bin/bash

#WITH INPUTFILE
#inputvalues.dat
FILE=$(realpath "$1")


if [ -z "$2" ]; then
	echo "No name supplied"
	NAME=$(awk -F= '$1=="NAME"{print $2;exit}' $FILE)
else
	NAME="$2"
fi

if [ -z "$3" ]; then
	echo "No name supplied"
	START=$(awk -F= '$1=="START"{print $2;exit}' $FILE)
else
	START=$(realpath "$3")
fi

CLUSTERSTARTIME=$(date +%s.%N)
printf "Cluster start: ${CLUSTERSTARTIME} for ${NAME}"

ROUND=$(awk -F= '$1=="ROUND"{print $2;exit}' $FILE)
ROUNDS=$(awk -F= '$1=="ROUNDS"{print $2;exit}'  $FILE)
PROG=$(awk -F= '$1=="PROGS"{print $2;exit}' $FILE)
PROGS="${PROG}/src"
SIMRNA=$(awk -F= '$1=="SIMRNA"{print $2;exit}'  $FILE)
WHERE=$(awk -F= '$1=="WHERE"{print $2;exit}' $FILE)
BUFFER=$(awk -F= '$1=="BUFFER"{print $2;exit}'  $FILE)
TYPE=$(awk -F= '$1=="TYPE"{print $2;exit}'  $FILE)
EXTEND=$(awk -F= '$1=="EXTEND"{print $2;exit}' $FILE)
RELAX=$(awk -F= '$1=="RELAX"{print $2;exit}' $FILE)
SIMROUND=$(awk -F= '$1=="SIMROUND"{print $2;exit}' $FILE)
SEED=$(awk -F= '$1=="SEED"{print $2;exit}' $FILE)
TREESEARCH=$(awk -F= '$1=="TREESEARCH"{print $2;exit}' $FILE)
CONTSEARCH1=$(awk -F= '$1=="CONTSEARCH1"{print $2;exit}' $FILE)
CONTSEARCH2=$(awk -F= '$1=="CONTSEARCH2"{print $2;exit}' $FILE)
CONSECUTIVEPERFECT=$(awk -F= '$1=="CONSECUTIVEPERFECT"{print $2;exit}' $FILE)
TARGET=$(awk -F= '$1=="TARGET"{print $2;exit}' $FILE)
EXPANDBMODE=$(awk -F= '$1=="EXPANDBMODE"{print $2;exit}'  $FILE)
SEQ="${START}/${NAME}.seq"

echo $NAME
echo $START
echo $ROUND
echo $ROUNDS
echo $PROGS
echo $SIMRNA
echo $WHERE
echo $BUFFER
echo $EXTEND
echo $RELAX
echo $SIMROUND


#####EXPANSION PREPERATION
#expansion to provide a bias through the SimRNA expansion
echo "EXPANTION"

# EXPANDBMODE 0: both directions at once
# EXPANDBMODE 1: only 'right'
# EXPANDBMODE 2: only 'left'
# EXPANDBMODE 3: alternate 'right' and 'left'
# EXPANDBMODE 4: alternate 'left' and 'right'
# EXPANDBMODE 5: first 'right' then 'left' ; 1 and then 2
# EXPANDBMODE 6: first 'left' then 'right' ; 2 and then 1
# EXPANDBMODE 7: user provided dotbracket notation with all intermediates

if [ "$EXPANDBMODE" -lt 7 ]; then
	for CURRENTROUND in `seq 1 1 ${ROUNDS}`; do

		ROLD="$(($CURRENTROUND-"1"))"
		NAMESS=$START/"${NAME}_${CURRENTROUND}.ss"
		NAMESSOLD=$START/"${NAME}_${ROLD}.ss"
		NAMESSTARGET=$START/"${NAME}_target.ss"

	  if [ "$TARGET" == true ] ; then
	    if [ "$EXPANDBMODE" = 0 ] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -l -o $NAMESS -t $NAMESSTARGET

	    elif [[ "$EXPANDBMODE" = 1 || "$EXPANDBMODE" = 5 ]] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS -t $NAMESSTARGET

	    elif [[ "$EXPANDBMODE" == 2 || "$EXPANDBMODE" == 6 ]] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS -t $NAMESSTARGET

	    elif [[ "$EXPANDBMODE" == 3 || "$EXPANDBMODE" == 4 ]]; then
	      if [ "$(( $CURRENTROUND % 2 ))" -eq 0 ]; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS -t $NAMESSTARGET
	      else
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS -t $NAMESSTARGET
	      fi

		  fi

		elif [ "$TARGET" == false ] ; then
	    if [ "$EXPANDBMODE" = 0 ] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -l -o $NAMESS

	    elif [[ "$EXPANDBMODE" = 1 || "$EXPANDBMODE" = 5 ]] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS

	    elif [[ "$EXPANDBMODE" = 2 || "$EXPANDBMODE" = 6 ]] ; then
	      python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS

	    elif [[ "$EXPANDBMODE" = 3 || "$EXPANDBMODE" = 4 ]]; then
	      if [ "$(( $CURRENTROUND % 2 ))" -eq 0 ]; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS
	      else
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS
	      fi

		  fi
	  fi

		if ([ ! -f $NAMESS ] && [ "$EXPANDBMODE" == 5 ]) || ([ ! -f $NAMESS ] && [ "$EXPANDBMODE" == 6 ]) ; then
	    if [ "$EXPANDBMODE" == 5 ] ; then
	      if [ "$TARGET" == true ] ; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS -t $NAMESSTARGET
	      elif [ "$TARGET" == false ] ; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -l -o $NAMESS
	      fi

	    elif [ "$EXPANDBMODE" == 6 ] ; then
	      if [ "$TARGET" == true ] ; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS -t $NAMESSTARGET
	      elif [ "$TARGET" == false ] ; then
	        python $PROGS/expandinteraction.py -b $BUFFER -n $SEQ -d $NAMESSOLD -r -o $NAMESS
	      fi

	    fi

	    if [ ! -f $NAMESS ] ; then
	      echo "New number of rounds: $ROLD"
	      ROUNDS="${ROLD}"
	      break 1
	    fi

	  elif [ ! -f $NAMESS ] ; then
			echo "New number of rounds: $ROLD"
			ROUNDS="${ROLD}"
			break 1
		fi

	done
fi

if ! [ ${RELAX} = "0" ]; then
	#####START from ernwin reconstructed pdb#####
	RELAXTIME=$(date +%s.%N)

	echo "$PROGS/SimRNA_config/job_simrna_start_${TYPE}.sh" $START $NAME $ROUND $SIMRNA $RELAX "${WHERE}" "${SEED}" $SIMROUND "$PROGS/SimRNA_config/"
	"$PROGS/SimRNA_config/job_simrna_start_${TYPE}.sh" $START $NAME $ROUND $SIMRNA $RELAX "${WHERE}" "${SEED}" $SIMROUND "$PROGS/SimRNA_config/"
	wait

	if grep -q "error(s)" "$START/${NAME}_${RELAX}_${ROUND}_1_1.log"; then
		cat "$START/${NAME}_${RELAX}_${ROUND}_1_1.log"
		echo "Error with ernwin PDB - want to change something?"
		read -n 1 -p Continue?
		"$PROGS/SimRNA_config/job_simrna_start_${TYPE}.sh" $START $NAME $ROUND $SIMRNA $RELAX "${WHERE}" "${SEED}" $SIMROUND "$PROGS/SimRNA_config/"
	fi

	mkdir $START/$ROUND

	mv "${START}/${NAME}_${RELAX}_${ROUND}"* ${START}/${ROUND}/.
	ROUNDNEW="$(($ROUND+"1"))"
	mkdir $START/${ROUNDNEW}
	cp "$START/${NAME}_0"* ${START}/${ROUND}/

	for step in $(seq 1 ${SIMROUND}); do
		TRAFL="${NAME}_${RELAX}_${ROUND}_${step}_${step}.trafl"
		PDB="${NAME}_${RELAX}_${ROUND}_${step}_${step}-000001.pdb"

		if [ "$TREESEARCH" == true ] ; then
			#write and analyse from every simRNA run all paralell runs (SimRounds)
			NAMESS=$START/$ROUND/"${NAME}_${ROUND}.ss"
			NAMECC=$START/$ROUND/"${NAME}_${ROUND}.ss_cc"
			NAMESTART=$START/"${NAME}_0.ss"

			if [ "$TARGET" == true ] ; then
				NAMEEND=$START/"${NAME}_target.ss"
			elif [ "$TARGET" == false ] ; then
				NAMEEND=$START/"${NAME}_${ROUNDS}.ss"
			fi

			mkdir ${START}/$ROUND/$step
			mkdir "${START}/$ROUND/$step/analyse"

			cp $START/$ROUND/${TRAFL} ${START}/$ROUND/${step}/analyse/
			cp ${START}/$ROUND/${PDB} ${START}/$ROUND/${step}/analyse/

			SimRNA_trafl2pdbs $START/$ROUND/${step}/analyse/*.pdb $START/$ROUND/${step}/analyse/*.trafl :
			wait
			python $PROGS/comparison.py -p $START/$ROUND/${step}/analyse/ -i $NAMECC -c $NAMESS -o ${NAME}_${RELAX}_0_${step}.csv -u ${NAME}_${RELAX}_0.csv -m 'w' -t $START/$ROUND/${step}/analyse/$TRAFL -s $NAMESTART -e $NAMEEND
			echo python $PROGS/comparison.py -p $START/$ROUND/${step}/analyse/ -i $NAMECC -c $NAMESS -o ${NAME}_${RELAX}_0_${step}.csv -u ${NAME}_${RELAX}_0.csv -m 'w' -t $START/$ROUND/${step}/analyse/$TRAFL -s $NAMESTART -e $NAMEEND
			wait
			mv ${NAME}_${RELAX}* ${START}/$ROUND/$step

			python ${PROGS}/selectnext.py -p ${START}/$ROUND/$step --print --first "${NAME}_${ROUND}.ss" --second "${NAME}_${ROUND}.ss_cc" -i "${NAME}_${RELAX}_${ROUND}" --${CONTSEARCH1} --consecutive $CONSECUTIVEPERFECT
			#get the forced structure and search within the individual runs
			wait
			SSCCBP=(${START}/"${NAME}_${RELAX}_0_${step}_${step}-******.ss_detected.bp") #e.g. CopStemsdesign1c0_expand_long03_0_2_2-000613.trafl
			BASIS=(${SSCCBP##*/})
			NR="$(cut -d'-' -f2 <<<${BASIS})"
			NR="$(cut -d'.' -f1 <<<${NR})"
			ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
			ROUNDSEL="$(cut -d'-' -f1 <<<${ROUNDSEL})"
			echo $BASIS
			echo $NR
			echo $ROUNDSEL

			SimRNA_trafl2pdbs "$START/$ROUND/$step/analyse/${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" "$START/$ROUND/$step/analyse/${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA
			wait
			mkdir $START/${ROUNDNEW}/${step} #for each of the simRNA-rounds

			cp $START/$ROUND/${step}/analyse/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}_AA.pdb" $START/$ROUNDNEW/${step}/"${NAME}_${ROUND}.pdb"
			cp $START/$ROUND/${step}/analyse/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected" $START/$ROUNDNEW/${step}/"${NAME}_${ROUND}.ss_cc"

			# interaction lenght + bp
			cp "$START/${NAME}_${ROUND}.il" $START/$ROUND/${step}/
			cp $START/$ROUND/"${NAME}_${ROUND}.il" ${START}/${ROUNDNEW}/${step}/"${NAME}_${ROUND}.il"
			mv "${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.il" ${START}/${ROUNDNEW}/${step}/"${NAME}_${ROUND}_ernwinrelax.il"
			mv "${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.bp" ${START}/${ROUNDNEW}/${step}/"${NAME}_${ROUND}.bp"

			#clean up directory
			rm -r "${START}/${ROUND}/${step}/analyse"

		elif [ "$TREESEARCH" == false ] ; then
			NAMESS=$START/"${NAME}_${ROUND}.ss"
			NAMECC=$START/"${NAME}_${ROUND}.ss_cc"
			NAMESTART=$START/"${NAME}_0.ss"

			if [ "$TARGET" == true ] ; then
				NAMEEND=$START/"${NAME}_target.ss"
			elif [ "$TARGET" == false ] ; then
				NAMEEND=$START/"${NAME}_${ROUNDS}.ss"
			fi
			
			mkdir $START/$ROUND/analyse

			mkdir ${START}/$ROUND/analyse/${step}

			cp $START/$ROUND/${TRAFL} ${START}/$ROUND/analyse/${step}/
			cp ${START}/$ROUND/${PDB} ${START}/$ROUND/analyse/${step}/

			SimRNA_trafl2pdbs $START/$ROUND/analyse/${step}/*.pdb $START/$ROUND/analyse/${step}/*.trafl :

			if [ $step = "1" ]; then
				python $PROGS/comparison.py -p $START/$ROUND/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${RELAX}_0_${step}.csv -u ${NAME}_${RELAX}_0.csv -m 'w' -t $START/$ROUND/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
			else #if no "treesearch" partI
				python $PROGS/comparison.py -p $START/$ROUND/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${RELAX}_0_${step}.csv -u ${NAME}_${RELAX}_0.csv -m 'a' -t $START/$ROUND/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
			fi

		else
			echo "Error wrong TREESEARCH value"
			exit 1
		fi
		wait
	done

	echo $NAME
	RELAXDURATION=$(echo "$(date +%s.%N) - $RELAXTIME" | bc)
	printf "Relaxationtime: %.6f seconds" $RELAXDURATION

	if [ "$TREESEARCH" == false ] ; then
		mv ${NAME}_${RELAX}* ${START}/$ROUND

		python ${PROGS}/selectnext.py -p ${START}/$ROUND --print --first "${NAME}_${ROUND}.ss" --second "${NAME}_${ROUND}.ss_cc" -i "${NAME}_${RELAX}_${ROUND}" --${CONTSEARCH1} --consecutive $CONSECUTIVEPERFECT
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

		SimRNA_trafl2pdbs $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA

		cp $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}_AA.pdb" $START/$ROUNDNEW/"${NAME}_${ROUND}.pdb"
		cp $START/$ROUND/analyse/$ROUNDSEL/"${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected" $START/$ROUNDNEW/"${NAME}_${ROUND}.ss_cc"

		# interaction lenght + bp
		mv $START/$ROUND/"${NAME}_${ROUND}.il" ${START}/${ROUNDNEW}/${step}/"${NAME}_${ROUND}.il"
		mv "${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.il" ${START}/${ROUNDNEW}/"${NAME}_${ROUND}_ernwinrelax.il"
		mv "${NAME}_${RELAX}_${ROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.bp" ${START}/${ROUNDNEW}/"${NAME}_${ROUND}.bp"

		#clean up directory
		cp "$START/${NAME}_${ROUND}.il" $START/$ROUND/
		mv "$START/${NAME}_${RELAX}"* ${START}/${ROUND}/
		rm -r "${START}/${ROUND}/analyse"
	fi #if no "treesearch" partII

	mv $START/$ROUND $START/"${ROUND}ernwinfinegrain"
	mv $START/$ROUNDNEW $START/${ROUND}

elif ([ $RELAX = "0" ] && [ "$TREESEARCH" == true ]) ; then #start without ernwin finegraining only for treesearch available
	mkdir $START/$ROUND
	ROUNDNEW="$(($ROUND+"1"))"
	mkdir $START/${ROUNDNEW}

	for simstep in $(seq 1 ${SIMROUND}); do #each of the starting SimRounds
		mkdir $START/$ROUND/${simstep}
		cp $START/"${NAME}_${ROUND}.il" $START/$ROUND/${simstep}/"${NAME}_${ROUND}.il"
		cp $START/"${NAME}_${ROUND}.ss_cc" $START/$ROUND/${simstep}/"${NAME}_${ROUND}.ss_cc"
		cp $START/"${NAME}_${ROUND}.pdb" $START/$ROUND/${simstep}/"${NAME}_${ROUND}.pdb"
	done
fi

#####EXPAND#####
for CURRENTROUND in `seq 0 1 ${ROUNDS}`; do
	ROUNDTIME=$(date +%s.%N)

	if [ "$TREESEARCH" == true ] ; then
		if [ "$CURRENTROUND" -ne 0 ]; then
			ROLD="$(($CURRENTROUND-"1"))"

			#if neither of the simsteps are expanded in the last round - stop simulation
			if [[ ! -z "$(ls -A $START/$CURRENTROUND)" ]]; then
    			echo "Directory is NOT empty - continue interaction expansion!"
			else
    			echo "Directory is empty! - no SimRNA-simulation can be expanded, end simulation!"
				rm $START/${CURRENTROUND}
				break 1
			fi
		fi

		for simstep in $(seq 1 ${SIMROUND}); do #each of the starting SimRounds
			echo "SimRound $simstep"

			if [ -f "$START/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.il" ]; then
				echo "File exist in Bash: "$START/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.il""

				echo "Compare only interaction"
				declare -i ssjet=$(cat "$START/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.il" ) #typeset -i
				declare -i ssold=$(cat "${START}/${ROLD}/${simstep}/${NAME}_${ROLD}.il")
				echo "Interaction length old: $ssold; Interaction length new: $ssjet"


				if ([ $CURRENTROUND -eq 0 ] && [ $ssjet -eq $ssold ]) || ([ $ssjet -gt $ssold ] && [ $CURRENTROUND -ne 0 ]) || ([ $CURRENTROUND -eq 0 ] && [ $RELAX = "0" ]); then #greater than -gt, equal -eq, not euqal -ne, || or, && and
					echo "Not the same = continue"

					cp "$START/${NAME}_${CURRENTROUND}.ss" "$START/${CURRENTROUND}/."

					cp "$START/${NAME}_${CURRENTROUND}.ss" "${START}/${CURRENTROUND}/${simstep}/."
					#expand SS
					NAMESS=${START}/${CURRENTROUND}/"${NAME}_${CURRENTROUND}.ss"
					NAMESSOLD=$START/${ROLD}/"${NAME}_${ROLD}.ss"

					#start the SimRNA expansion job
					"$PROGS/SimRNA_config/job_simrna_start_${TYPE}.sh" ${START}/${CURRENTROUND}/${simstep}/ $NAME $CURRENTROUND $SIMRNA $EXTEND $WHERE $SEED $SIMROUND "$PROGS/SimRNA_config/"
					wait

					mkdir ${START}/${CURRENTROUND}/${simstep}/analyse

					for step in $(seq 1 ${SIMROUND}); do
						mkdir ${START}/${CURRENTROUND}/${simstep}/analyse/${step}
						TRAFL="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}.trafl"
						PDB="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}-000001.pdb"
						cp ${START}/${CURRENTROUND}/${simstep}/${TRAFL} ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/
						cp ${START}/${CURRENTROUND}/${simstep}/${PDB} ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/
						SimRNA_trafl2pdbs ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/*.pdb ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/*.trafl :
						wait

						NAMESTART=$START/"${NAME}_0.ss"
						if [ "$TARGET" == true ] ; then
							NAMEEND=$START/"${NAME}_target.ss"
						elif [ "$TARGET" == false ] ; then
							NAMEEND=$START/"${NAME}_${ROUNDS}.ss"
						fi

						if [ $step = "1" ]; then
							#first round + create outputfiles
							python $PROGS/comparison.py -p ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'w' -t ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
							echo python $PROGS/comparison.py -p ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'w' -t ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
							wait
						else
							#append the following output
							python $PROGS/comparison.py -p ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'a' -t ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
							echo python $PROGS/comparison.py -p ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/ -i ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss_cc -c ${START}/${CURRENTROUND}/${simstep}/${NAME}_${CURRENTROUND}.ss -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'a' -t ${START}/${CURRENTROUND}/${simstep}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
							wait
						fi
						wait
					done

					mv ${NAME}_${EXTEND}* $START/${CURRENTROUND}/${simstep}/

					#find the starting structure for the next round
					python ${PROGS}/selectnext.py -p ${START}/${CURRENTROUND}/${simstep} --print --first "${NAME}_${CURRENTROUND}.ss" --second "${NAME}_${CURRENTROUND}.ss_cc" --${CONTSEARCH2} -i "${NAME}_${EXTEND}_${CURRENTROUND}" --consecutive $CONSECUTIVEPERFECT
					wait
					#BASIS=$(basename -s .ss_detected.bp "${START}/${CURRENTROUND}/${simstep}/*.ss_detected.bp")
					SSCCBP=($START/${CURRENTROUND}/${simstep}/*.ss_detected.bp)
					BASIS=(${SSCCBP##*/})
					NR="$(cut -d'-' -f2 <<<${BASIS})"
					NR="$(cut -d'.' -f1 <<<${NR})"
					ROUNDSEL="$(cut -d'_' -f5 <<<${BASIS})"
					ROUNDSEL="$(cut -d'-' -f1 <<<${ROUNDSEL})"
					echo $BASIS
					echo $NR
					echo $ROUNDSEL

					#calculate a full atom strucutre
					SimRNA_trafl2pdbs $START/$CURRENTROUND/${simstep}/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-000001.pdb" $START/$CURRENTROUND/${simstep}/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}.trafl" ${NR} AA
					wait
					#prepare everything for the next round
					ROUNDNEW="$(($CURRENTROUND+"1"))"
					mkdir $START/$ROUNDNEW
					mkdir $START/$ROUNDNEW/${simstep}

					NEWSSCC=$START/${CURRENTROUND}/${simstep}/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected"
					cp $NEWSSCC ${START}/${ROUNDNEW}/${simstep}/"${NAME}_${ROUNDNEW}.ss_cc"

					NEWPDB=$START/$CURRENTROUND/${simstep}/analyse/$ROUNDSEL/"${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}_AA.pdb"
					mv $NEWPDB $START/$ROUNDNEW/${simstep}/"${NAME}_${ROUNDNEW}.pdb"

					# interaction lenght + bp
					mv "${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.il" ${START}/${ROUNDNEW}/${simstep}/"${NAME}_${ROUNDNEW}.il"
					mv "${NAME}_${EXTEND}_${CURRENTROUND}_${ROUNDSEL}_${ROUNDSEL}-${NR}.ss_detected.bp" ${START}/${ROUNDNEW}/${simstep}/"${NAME}_${ROUNDNEW}.bp"

					#clean up directory
					rm -r "$START/${CURRENTROUND}/${simstep}/analyse"

				else
					echo "Interaction length the same stop here"
				fi

			else
				echo "Interaction il-files does not exist"
			fi

			wait
		done


	elif [ "$TREESEARCH" == false ] ; then
		ROLD="$(($CURRENTROUND-"1"))"
		cp "$START/${NAME}_${CURRENTROUND}.ss" "$START/${CURRENTROUND}/."
		#expand SS
		NAMESS=$START/${CURRENTROUND}/"${NAME}_${CURRENTROUND}.ss"
		NAMESSOLD=$START/${ROLD}/"${NAME}_${ROLD}.ss"
		NAMECC=${START}/${CURRENTROUND}/"${NAME}_${CURRENTROUND}.ss_cc"
		NAMESTART=$START/"${NAME}_0.ss"

		if [ "$TARGET" == true ] ; then
	    	NAMEEND=$START/"${NAME}_target.ss"
	    elif [ "$TARGET" == false ] ; then
	        NAMEEND=$START/"${NAME}_${ROUNDS}.ss"
	    fi

		#+ checking algorithm - now both  before the SimRNA expantion

		#start the SimRNA expansion job
		"$PROGS/SimRNA_config/job_simrna_start_${TYPE}.sh" $START/${CURRENTROUND}/ $NAME $CURRENTROUND $SIMRNA $EXTEND $WHERE $SEED $SIMROUND "$PROGS/SimRNA_config/"
		wait

		#analyse of the SimRNA run
		mkdir $START/${CURRENTROUND}/analyse
		for step in $(seq 1 ${SIMROUND}); do
			mkdir ${START}/${CURRENTROUND}/analyse/${step}
			TRAFL="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}.trafl"
			PDB="${NAME}_${EXTEND}_${CURRENTROUND}_${step}_${step}-000001.pdb"
			cp $START/${CURRENTROUND}/${TRAFL} ${START}/${CURRENTROUND}/analyse/${step}/
			cp ${START}/${CURRENTROUND}/${PDB} ${START}/${CURRENTROUND}/analyse/${step}/
			SimRNA_trafl2pdbs $START/${CURRENTROUND}/analyse/${step}/*.pdb $START/${CURRENTROUND}/analyse/${step}/*.trafl :

			if [ $step = "1" ]; then
				#first round + create outputfiles
				python $PROGS/comparison.py -p $START/${CURRENTROUND}/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'w' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
			else
				#append the following output
				python $PROGS/comparison.py -p $START/${CURRENTROUND}/analyse/${step}/ -i $NAMECC -c $NAMESS -o ${NAME}_${EXTEND}_${CURRENTROUND}_${step}.csv -u ${NAME}_${EXTEND}_${CURRENTROUND}.csv -m 'a' -t $START/${CURRENTROUND}/analyse/${step}/$TRAFL -s $NAMESTART -e $NAMEEND
			fi
		done

		mv ${NAME}_${EXTEND}* ${START}/${CURRENTROUND}/

		#find the starting structure for the next round
		python ${PROGS}/selectnext.py -p ${START}/${CURRENTROUND} --print --first "${NAME}_${CURRENTROUND}.ss" --second "${NAME}_${CURRENTROUND}.ss_cc" --${CONTSEARCH2} -i "${NAME}_${EXTEND}_${CURRENTROUND}" --consecutive $CONSECUTIVEPERFECT

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
		rm -r "${START}/${CURRENTROUND}/analyse"


		echo "Compare only interaction"
		declare -i ssjet=$(cat "${START}/${CURRENTROUND}/${NAME}_${CURRENTROUND}.il" ) #typeset -i
		declare -i ssnew=$(cat "${START}/${ROUNDNEW}/${NAME}_${ROUNDNEW}.il")
		echo "Interaction length old: $ssold; Interaction length new: $ssjet"

		if [[ $ssnew -gt $ssjet ]]; then #greater than -gt
			echo "Not the same = continue"
		else
			echo "The same"
			#rm -r data
			exit 0
		fi

	else
		echo "error with the expansion loop"
		exit 1

	fi

	echo $NAME $CURRENTROUND
	ROUNDDURATION=$(echo "$(date +%s.%N) - $ROUNDTIME" | bc)
	printf "Execution time for this extension: %.6f seconds" $ROUNDDURATION
done

DURATIONCLUSTER=$(echo "$(date +%s.%N) - $CLUSTERSTARTIME" | bc)
printf "Execution time for this cluster: %.6f seconds" $DURATIONCLUSTER

echo "SIMULATION DONE"

#rm -r data
