#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

decomposePar

mpirun -np 2 lammpsFoamParallelDebug -parallel &> log.parallel

reconstructPar -latestTime
