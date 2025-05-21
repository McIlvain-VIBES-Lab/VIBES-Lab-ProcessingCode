#!/bin/bash
export FSLDIR=/Users/vibes/fsl
source $FSLDIR/etc/fslconf/fsl.sh

fslorient -setqformcode 1 "$1"
fslorient -copyqform2sform "$1"
fslmodhd "$1" scl_slope 1