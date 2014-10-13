#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

blockMesh
decomposePar

mpirun -np 4 lammpsFoam -parallel &> log.parallel

reconstructPar -latestTime
