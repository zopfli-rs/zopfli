#!/bin/bash
echo 'SHA,Lines of C,Lines of Rust' > results.csv
START=$(git rev-parse --abbrev-ref HEAD)
for COMMIT in $(git rev-list upstream/master.. --reverse)
do
    git checkout $COMMIT

    git ls-files --error-unmatch src/zopfli/*.{c,h} > /dev/null 2>&1
    C_FILES_PRESENT=$?
    if [ $C_FILES_PRESENT -eq 0 ]
    then
        C=$(wc -l $(git ls-files src/zopfli/*.{c,h}) | tail -n 1 | awk '{print $1}')
    else
        C=0
    fi

    git ls-files --error-unmatch *.rs > /dev/null 2>&1
    RUST_FILES_PRESENT=$?
    if [ $RUST_FILES_PRESENT -eq 0 ]
    then
        RUST=$(wc -l $(git ls-files *.rs) | tail -n 1 | awk '{print $1}')
    else
        RUST=0
    fi

    echo "$COMMIT,$C,$RUST" >> results.csv
done
git checkout $START
