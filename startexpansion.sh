#!/bin/bash

## ./startexpansion.sh inputvalues.dat

FILE=$(realpath "$1")
CLUSTER=$(awk -F= '$1=="CLUSTER"{print $2;exit}' $FILE)
NAME=$(awk -F= '$1=="NAME"{print $2;exit}' $FILE)
START=$(awk -F= '$1=="START"{print $2;exit}' $FILE)


if [ $CLUSTER = "0" ]; then
  ./main.sh $FILE "${NAME}" "$START"

else
  CLUSTER="$(($CLUSTER-"1"))"

  for CST in `seq 0 1 ${CLUSTER}`; do
    cp "$START/${NAME}_0.il" "$START/cluster${CST}/${NAME}c${CST}_0.il" #check in startscript!
    ./main.sh $FILE "${NAME}c${CST}" "$START/cluster${CST}" 
    wait
  done

fi
