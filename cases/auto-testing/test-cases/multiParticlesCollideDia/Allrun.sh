#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

f21x

blockMesh
lammpsFoam

./particlePosition.py
./particleVelocity.py
