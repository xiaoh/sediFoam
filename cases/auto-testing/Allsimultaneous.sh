#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDIR=$PWD
reportDIR=$PWD/test-report-generation


codeVersion=`git log -1 --format="%h"`
cd $reportDir

cat body/before.tex > lammps-paper.tex
reportName="report"
reportName+=$codeVersion

mkdir $reportDIR/$reportName
cd $reportDIR/$reportName
mkdir figs

generationDIR=$PWD
cd $reportDIR/report-example
cp -rf * $generationDIR

command="pwd"
index=1
cd $currentDIR/test-cases
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

cd $currentDIR/test-cases
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    cd $generationDIR/figs
    cp -rf $xDIR/data/*.pdf . 
    cd $xDIR
    cd $benchDIR
done

cd $generationDIR
./generate.sh
