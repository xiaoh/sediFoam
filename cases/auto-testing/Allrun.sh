#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

codeVersion="ab04"
currentDIR=$PWD
reportDIR=$PWD/reports

cd $currentDIR/bench
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    ./Allrun.sh &> log.run
    if [ ! -d $reportDIR/figs ]; then
	mkdir -p $reportDIR/figs
    fi
    cd $reportDIR/figs
    cp -rf $xDIR/data/*.pdf . 
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
    cp -rf $xDIR/data/*.pdf . 
    cd $xDIR
    cd $exampleDIR
done

cd $reportDIR
./generate.sh
