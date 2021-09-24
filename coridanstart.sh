#!/bin/bash -x

## ./coridanstart.sh inputvalues_100-000.dat

FILE=$(realpath "$1")
CLUSTER=$(awk -F= '$1=="CLUSTER"{print $2;exit}' $FILE)
NAME=$(awk -F= '$1=="NAME"{print $2;exit}' $FILE)
START=$(awk -F= '$1=="START"{print $2;exit}' $FILE)


for CST in `seq 0 1 ${CLUSTER}`; do
  ./ultimatescript.sh $FILE "${NAME}c${CST}" "$START/cluster${CST}"
done
