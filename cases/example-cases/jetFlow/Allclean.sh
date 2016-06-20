#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

rm -rf 0/polyMesh
rm -rf processor*
rm -rf *[1-9]*

rm -rf log*
rm -rf snap*
