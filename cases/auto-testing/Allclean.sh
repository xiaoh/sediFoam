#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDIR=$PWD
reportDIR=$PWD/test-report-generation

cd $currentDIR/test-cases
benchDIR=$PWD
for folder in *; do
    cd $folder 
    echo "cleaning: $PWD"
    ./Allclean.sh 
    rm data/*.pdf
    cd $benchDIR
done

cd $currentDIR/test-cases
exampleDIR=$PWD
for folder in *; do
    cd $folder 
    echo "cleaning: $PWD"
    ./Allclean.sh 
    rm data/*.pdf
    cd $exampleDIR
done

# cd $reportDIR
# 
# for folder in *; do
#     cd $folder 
#     echo "cleaning: $PWD"
#     rm *.log *.pdf *.aux *.out
#     cd ..
# done
