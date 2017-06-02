#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

blockMesh > log.blockMesh
decomposePar > log.decomposePar
mpirun -np 48 lammpsFoam -parallel > log.parallel
