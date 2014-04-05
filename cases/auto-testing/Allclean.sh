#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDIR=$PWD
reportDIR=$PWD/report

cd $currentDIR/bench
benchDIR=$PWD
for folder in *; do
    cd $folder 
    echo "cleaning: $PWD"
    ./Allclean.sh 
    rm data/*.pdf
    cd $benchDIR
done

cd $currentDIR/example
exampleDIR=$PWD
for folder in *; do
    cd $folder 
    echo "cleaning: $PWD"
    ./Allclean.sh 
    rm data/*.pdf
    cd $exampleDIR
done

cd $reportDIR
rm figs/*.pdf
