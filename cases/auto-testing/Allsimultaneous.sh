#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Save directory names
currentDIR=$PWD
reportDIR=$PWD/test-report-generation

# Set up a new DIR with the name of the code version
codeVersion=`git log -1 --format="%h"`
cd $reportDir

reportName="report-simu-"
reportName+=$codeVersion

mkdir $reportDIR/$reportName
cd $reportDIR/$reportName
mkdir figs

# Copy necessary files to generate the report
generationDIR=$PWD
cd $reportDIR/essential
cp -rf * $generationDIR

# Run all cases simultaneously
command="pwd"
index=1
cd $currentDIR/test-cases
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    cd $benchDIR
    index=$[index+1]
    command+=" & sh "
    command+=$xDIR/Allrun.sh
done

cd $currentDIR
eval $command
wait

# Copy necessary figures to generate the report
cd $currentDIR/test-cases
benchDIR=$PWD
for folder in *; do
    cd $folder 
    xDIR=$PWD
    cd $generationDIR/figs
    cp -rf $xDIR/data/*.pdf . 
    cd $benchDIR
done

# Generate the report
cd $generationDIR
./generate.sh
