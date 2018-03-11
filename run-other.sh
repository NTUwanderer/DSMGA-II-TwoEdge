#!/bin/bash
REPEAT=100
nThreads=$((`nproc --all` - 2))
if [ "$#" == "0" ]; then
echo "./run-other nInitial FFF DIR ell1 ell2 ... ellN"
exit 1
fi
nInitial=$1
FFF=$2
DIR=$3
mkdir $DIR
while [ $# -gt  3 ]
do
    for (( i=0; $i<$REPEAT; i=$i+1 ))
    do
        NUM=`echo $i | awk '{printf "%03d", $1}'`
        ELL=$4
        nohup nice -n19 ./DSMGA2 $ELL $nInitial $FFF 10000 -1 1 1 -1 > ./$DIR/$ELL-$NUM &
        echo "Submitting $ELL-$NUM"
        sleep 2
        TT=$(ps -xaf | grep DSMGA2 | wc -l)
        while [ $TT -gt $nThreads ]
        do
            sleep 1
            TT=$(ps -xaf | grep DSMGA2 | wc -l)
        done
    done
    shift
done