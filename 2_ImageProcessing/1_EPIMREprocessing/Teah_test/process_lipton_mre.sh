#!/bin/bash

#Script to make NIFTIs into the correct format for Lipton Lab 
#Note: We are now using this script for other studies 

export FSLDIR=/Users/vibes_1/fsl
source $FSLDIR/etc/fslconf/fsl.sh

fslorient -setqformcode 1 "$1"
fslorient -copyqform2sform "$1"
fslmodhd "$1" scl_slope 1