#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

./Allclean.sh

blockMesh
decomposePar

mpirun -np 3 lammpsFoam -parallel &> log.parallel

reconstructPar -latestTime

./particlePosition.py
./particleVelocity.py
