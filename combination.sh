#!/bin/bash -x

FILE=$(realpath "$1")
START=$(awk -F= '$1=="START"{print $2;exit}' $FILE)
DESIGNS=$(awk -F= '$1=="DESIGNS"{print $2;exit}' $FILE)
NAME=$(awk -F= '$1=="NAME"{print $2;exit}' $FILE)
BASENAME=$(awk -F= '$1=="BASENAME"{print $2;exit}' $FILE)
COARSE=$(awk -F= '$1=="ERNWIN"{print $2;exit}' $FILE)
PROGS=$(awk -F= '$1=="PROGS"{print $2;exit}' $FILE)
ERNROUND=$(awk -F= '$1=="ERNROUND"{print $2;exit}' $FILE)
ERNITERATIONS=$(awk -F= '$1=="ERNITERATIONS"{print $2;exit}' $FILE)
CLUSTER=$(awk -F= '$1=="CLUSTER"{print $2;exit}' $FILE)
FALLBACKSTATES=$(awk -F= '$1=="FALLBACKSTATES"{print $2;exit}' $FILE)

ERNWIN="$COARSE/fess/scripts/ernwin.py"
RECONSTRUCTION="$COARSE/fess/scripts/reconstruct.py"
CGS=/"$COARSE/RESOURCES/CGS"
PDB="$COARSE/RESOURCES/PDB_DIR/"

if [ "$FALLBACKSTATES" == true ] ; then
  FALLBACK="$COARSE/RESOURCES/fallback.stats"
fi


####SECONDDESIGN=(1 3 7 8 11 13 14 16 17 19 20 21 22 23 24 25 26 29 34 35 36 37 39 41 44 45 46 49 52 53 54 55 56 57 58 59 61 62 63 68 69 71 72 76 78 79 80 83 87 89 90 92 93 98 99 100)

######./ultimatescriptstart.sh "${START}/$DESIGN" "CopStemsdesign${DESIGN}" "/home/irene/Programs/GitHub/RNA-Interaction-Workflow/inputvalues.dat"

if [[ ${DESIGNS} == 0 ]]; then
  if [ "$CLUSTER" -eq "0" ]; then

    if [ "$FALLBACKSTATES" == true ] ; then
      python2 $ERNWIN "$START/${NAME}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots --fallback-stats $FALLBACK
    else
      python2 $ERNWIN "$START/${NAME}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots
    fi
    wait
    mv ${NAME} "$START/"
    COUNTER=0

    while [  $COUNTER -lt $ERNROUND ]; do
      echo "ernwin reconstruction $COUNTER"

      if [ "$FALLBACKSTATES" == true ] ; then
        python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/best$COUNTER.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
      else
        python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/best$COUNTER.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
      fi

      wait
      RECONSTRUCTIONFILE="$START/${NAME}/simulation_01/best$COUNTER.coord.reconstr.pdb"

      #check if file exist
      if [ ! -f "$RECONSTRUCTIONFILE" ]; then
        let COUNTER=$COUNTER+1
      else
        echo "possible ernwin reconstruction $COUNTER"
        cp $RECONSTRUCTIONFILE  "$START/${NAME}_0.pdb"
        ./ultimatescript.sh $FILE "${NAME}" "$START"

        break 1
      fi

      echo "ERROR: No ernwin reconstruction possible"

      retVal=$?
      TRY="0"
      ERNWINCOUNT="1"

      while [[ $retVal -ne 0 && "$TRY" -lt $ERNROUND ]]; do
        echo "Error in PDB translation"

        while [ $ERNWINCOUNT -lt $ERNROUND ]; do
          #ERNWINCOUNT="$(($COUNTER+$RECONCOUNT))"
          if [ "$FALLBACKSTATES" == true ] ; then
            python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/best$ERNWINCOUNT.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
          else
            python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/best$ERNWINCOUNT.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
          fi
          wait
          RECONSTRUCTIONFILE="$START/${NAME}/simulation_01/best$ERNWINCOUNT.coord.reconstr.pdb"

          #check if file exist
          if [ ! -f "$RECONSTRUCTIONFILE" ]; then
            let COUNTER=$COUNTER+1
          else
            echo "possible ernwin reconstruction $COUNTER"
            cp $RECONSTRUCTIONFILE  "$START/${NAME}_0.pdb"
            break 1
          fi

          echo "ERROR: No ernwin reconstruction possible"
        done

        ./ultimatescript.sh $FILE "${NAME}" "$START"
        retVal=$?

        if [ $retVal -ne 0 ]; then
            echo "Error in PDB translation"
            let TRY=$TRY+1
            let ERNWINCOUNT=$ERNWINCOUNT+1
        fi

      done
    done

  #Calculation with ernwin clustering
  else
    if [ "$FALLBACKSTATES" == true ] ; then
      python2 $ERNWIN "$START/${NAME}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots --fallback-stats $FALLBACK
    else
      python2 $ERNWIN "$START/${NAME}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots
    fi
    wait
    mv ${NAME} "$START/"
    COUNTER=0

    #python $PROGS/"ernwindiversity.py" -i "$START/${NAME}/simulation_01/" -n $ERNROUND -c $CLUSTER

    for CST in `seq 0 1 ${CLUSTER}`; do
      while IFS= read -r line; do
        echo "From cluster $CST sample : $line"

        if [ "$FALLBACKSTATES" == true ] ; then
          python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/$line.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
        else
          python2 $RECONSTRUCTION "$START/${NAME}/simulation_01/$line.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
        fi

        mv ${NAME} "$START/"
        wait
        RECONSTRUCTIONFILE="$START/${NAME}/simulation_01/$line.coord.reconstr.pdb"

        #check if file exist
        if [ -f "$RECONSTRUCTIONFILE" ]; then
          echo "possible ernwin reconstruction $line"
          mkdir "$START/cluster${CST}"
          cp $RECONSTRUCTIONFILE "$START/cluster${CST}/${NAME}c${CST}_0.pdb"
          cp "$START/${NAME}.seq" "$START/cluster${CST}/${NAME}c${CST}.seq"
          cp "$START/${NAME}_0.ss" "$START/cluster${CST}/${NAME}c${CST}_0.ss"
          cp "$START/${NAME}_0.ss_cc" "$START/cluster${CST}/${NAME}c${CST}_0.ss_cc"
          cp "$START/${NAME}_0.il" "$START/cluster${CST}/${NAME}c${CST}_0.il"
          cp "$START/${NAME}_target.ss" "$START/cluster${CST}/${NAME}c${CST}_target.ss"

          ./ultimatescript.sh $FILE "${NAME}c${CST}" "$START/cluster${CST}"
          break 1
        fi

        echo "ERROR: No ernwin reconstruction possible"
      done < "$START/${NAME}/simulation_01/cluster$CST.csv"
    done
  fi

else
  python $PROGS/"RNAdesign.py" -i "$START/${BASENAME}_0.ss" -o ${BASENAME} -n $DESIGNS -s $DESIGNS #Scenario testing
  mv "${NAME}"* $START/.

  python "$PROGS/formattranslation.py" -p $START -n ${NAME} -c $DESIGNS

  #calculation without clustering of the ernwin structures
  if [ "$CLUSTER" -eq "0" ]; then

    for DESIGN in `seq 1 1 ${DESIGNS}`; do

      mkdir "$START/${NAME}${DESIGN}"
      cp "$START/${NAME}design${DESIGN}.seq" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.seq"
      cp "$START/${NAME}_0.ss" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss"
      cp "$START/${NAME}_0.ss_cc" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss_cc"
      cp "$START/${NAME}design${DESIGN}.fa" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa"
      cp "$START/${NAME}_0.il"  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.il"
      cp "$START/${NAME}_target.ss"  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_target.ss"

      if [ "$FALLBACKSTATES" == true ] ; then
        python2 $ERNWIN "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots --fallback-stats $FALLBACK
      else
        python2 $ERNWIN "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots
      fi

      wait
      mv ${NAME}${DESIGN} "$START/${NAME}${DESIGN}/"

      COUNTER=0

      while [  $COUNTER -lt $ERNROUND ]; do #Scenario testing
        echo "ernwin reconstruction $COUNTER"

        if [ "$FALLBACKSTATES" == true ] ; then
          python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$COUNTER.coord"  --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
        else
          python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$COUNTER.coord"  --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
        fi

        wait
        RECONSTRUCTIONFILE="$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$COUNTER.coord.reconstr.pdb"

        #check if file exist
        if [ ! -f "$RECONSTRUCTIONFILE" ]; then
          let COUNTER=$COUNTER+1
        else
          echo "possible ernwin reconstruction $COUNTER"
          cp $RECONSTRUCTIONFILE  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.pdb"
          break 1
        fi

        echo "ERROR: No ernwin reconstruction possible"
      done

      ./ultimatescript.sh $FILE "${NAME}design${DESIGN}" "$START/${NAME}${DESIGN}/"

      retVal=$?
      TRY="0"
      ERNWINCOUNT="1"

      while [[ $retVal -ne 0 && "$TRY" -lt $ERNROUND ]]; do
        echo "Error in PDB translation"

        while [ $ERNWINCOUNT -lt $ERNROUND ]; do
          #ERNWINCOUNT="$(($COUNTER+$RECONCOUNT))"
          if [ "$FALLBACKSTATES" == true ] ; then
            python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$ERNWINCOUNT.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
          else
            python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$ERNWINCOUNT.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
          fi

          wait
          RECONSTRUCTIONFILE="$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/best$ERNWINCOUNT.coord.reconstr.pdb"

          #check if file exist
          if [ ! -f "$RECONSTRUCTIONFILE" ]; then
            let COUNTER=$COUNTER+1
          else
            echo "possible ernwin reconstruction $COUNTER"
            cp $RECONSTRUCTIONFILE  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.pdb"
            break 1
          fi

          echo "ERROR: No ernwin reconstruction possible"
        done

        ./ultimatescript.sh $FILE "${NAME}design${DESIGN}" "$START/${NAME}${DESIGN}/"
        retVal=$?

        if [ $retVal -ne 0 ]; then
            echo "Error in PDB translation"
            let TRY=$TRY+1
            let ERNWINCOUNT=$ERNWINCOUNT+1
        fi

      done
    done

  #Calculation with ernwin clustering
  else
    for DESIGN in `seq 1 1 ${DESIGNS}`; do

      mkdir "$START/${NAME}${DESIGN}"
      cp "$START/${NAME}design${DESIGN}.seq" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.seq"
      cp "$START/${NAME}_0.ss" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss"
      cp "$START/${NAME}_0.ss_cc" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss_cc"
      cp "$START/${NAME}design${DESIGN}.fa" "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa"
      cp "$START/${NAME}_0.il"  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.il"
      cp "$START/${NAME}_target.ss"  "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_target.ss"

      if [ "$FALLBACKSTATES" == true ] ; then
        python2 $ERNWIN "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots --fallback-stats $FALLBACK
      else
        python2 $ERNWIN "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.fa" --save-n-best $ERNROUND --iterations $ERNITERATIONS --pseudoknots
      fi

      wait
      mv ${NAME}${DESIGN} "$START/${NAME}${DESIGN}/"

      COUNTER=0

      python $PROGS/"ernwindiversity.py" -i "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/" -n $ERNROUND -c $CLUSTER

      CLUSTER="$(($CLUSTER-"1"))"

      for CST in `seq 0 1 ${CLUSTER}`; do
        while IFS= read -r line; do
          echo "From cluster $CST sample : $line"
          if [ "$FALLBACKSTATES" == true ] ; then
            python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/$line.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken --fallback-stats $FALLBACK
          else
            python2 $RECONSTRUCTION "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/$line.coord" --source-cg-dir $CGS --source-pdb-dir $PDB --reassign-broken
          fi

          wait
          RECONSTRUCTIONFILE="$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/$line.coord.reconstr.pdb"

          #check if file exist
          if [ -f "$RECONSTRUCTIONFILE" ]; then
            echo "possible ernwin reconstruction $line"
            mkdir "$START/${NAME}${DESIGN}/cluster${CST}"
            cp $RECONSTRUCTIONFILE  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}_0.pdb"
            cp "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}.seq"  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}.seq"
            cp "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss"  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}_0.ss"
            cp "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.ss_cc"  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}_0.ss_cc"
            cp "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_0.il"  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}_0.il"
            cp "$START/${NAME}${DESIGN}/${NAME}design${DESIGN}_target.ss"  "$START/${NAME}${DESIGN}/cluster${CST}/${NAME}design${DESIGN}c${CST}_target.ss"

            ./ultimatescript.sh $FILE "${NAME}design${DESIGN}c${CST}" "$START/${NAME}${DESIGN}/cluster${CST}"
            break 1
          fi

          echo "ERROR: No ernwin reconstruction possible"
        done < "$START/${NAME}${DESIGN}/${NAME}${DESIGN}/simulation_01/cluster$CST.csv"
      done
    done
  fi
fi
