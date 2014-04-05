#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

rm -rf processor*
rm -rf *[1-9]*
rm -rf probes 
rm -rf sets
rm log.*
rm snapshot.*
