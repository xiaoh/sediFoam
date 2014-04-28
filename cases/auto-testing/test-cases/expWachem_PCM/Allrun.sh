#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

f21x

./AllrunParallel.sh

./postprocessing.py
