#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing lammpsFoam.."
currentDir=$PWD
echo -n "Enter the directory of your LAMMPS and press [ENTER]: "
read lammpsDir

# Determine if the directory of LAMMPS exists or not.
# If not, look for LAMMPS in the default directory.
if [ ! -d "$lammpsDir" ]
then
    echo "Directory NOT found! Use default directory instead."
    lammpsDir="$PWD/lammps-1Feb14"
fi

cd $lammpsDir
lammpsDir=$PWD

echo "Directory of LAMMPS is: " $lammpsDir

cd $lammpsDir/src
# Make packages
make clean-all
make no-all

cd $currentDir/lammpsFoam
rm Make/options

wclean dragModels
wclean 

