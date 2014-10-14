#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

cd constant/polyMesh
rm -rf boundary points faces neighbour owner
cd ../..

rm *[1-9]* -rf
rm log*
rm snap*
