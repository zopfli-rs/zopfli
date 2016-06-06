#!/bin/bash
echo 'SHA,result' > results.csv
START=$(git rev-parse --abbrev-ref HEAD)
for COMMIT in $(git rev-list upstream/master.. --reverse)
do
    git checkout $COMMIT
    make clean && make zopfli
    if [ $? -ne 0 ]
    then
        break
    fi
    # time=$( TIMEFORMAT="%R"; { time ./test/run.sh; } 2>&1 )
    # git checkout test/results
    # echo "$COMMIT,$time" >> results.csv
done
