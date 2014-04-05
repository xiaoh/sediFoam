#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDIR=$PWD
reportDIR=$PWD/report

command="pwd"
index=1
cd $currentDIR/bench
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    cd $benchDIR
    index=$[index+1]
    echo $index
    echo $xDIR
    command+=" & sh "
    command+=$xDIR/Allrun.sh
done

cd $currentDIR
eval $command
wait

cd $currentDIR/bench
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    cd $reportDIR/figs
    ln -sf $xDIR/data/*.pdf . 
    cd $xDIR
    cd $benchDIR
done

cd $currentDIR/example
exampleDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    ./Allrun.sh &> log.run
    cd $reportDIR/figs
    ln -sf $xDIR/data/*.pdf . 
    cd $xDIR
    cd $exampleDIR
done

cd $reportDIR
./generate.sh
