#!/bin/bash
echo 'SHA,result' > results.csv
START=$(git rev-parse --abbrev-ref HEAD)
for COMMIT in $(git rev-list upstream/master.. --reverse)
do
    git checkout $COMMIT
    if [ $? -ne 0 ]
    then
        echo "COULD NOT CHECKOUT"
        break
    fi
    make clean && make zopfli
    if [ $? -ne 0 ]
    then
        echo "COULD NOT MAKE"
        break
    fi
    time=$( TIMEFORMAT="%R"; { time ./test/run.sh; } 2>&1 )
    echo "$COMMIT,$time" >> results.csv
done
git checkout $START
