#!/bin/bash
REPEAT=100
nThreads=$((`nproc --all` - 2))
if [ "$#" == "0" ]; then
echo "./run-spin nInitial DIR ell1 ell2 ... ellN"
exit 1
fi
nInitial=$1
DIR=$2
mkdir -p $DIR
while [ $# -gt  2 ]
do
    for (( i=1; $i<=$REPEAT; i=$i+1 ))
    do
        NUM=`echo $i | awk '{printf "%03d", $1}'`
        ELL=$3
        nohup nice -n19 ./DSMGA2 $ELL $nInitial 5 10000 10000000 1 1 -1 $i > ./$DIR/$ELL-$NUM &
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
