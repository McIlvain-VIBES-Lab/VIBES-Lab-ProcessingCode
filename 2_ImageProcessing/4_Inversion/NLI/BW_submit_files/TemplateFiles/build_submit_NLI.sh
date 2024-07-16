#!/bin/bash
#
# Build runfile and PBS script, then submit PBS for NLI reconstruction
#
# ------------------------------------------------
#
subj=$1
#
nod=`ls $PWD/*.hom.nod`
elm=`ls $PWD/*.hom.elm`
bnd=`ls $PWD/*.bnd`
reg=`ls $PWD/*.regstack`
dsp=`ls $PWD/*.dsp`
invDir=$2
out="$PWD/$invDir/MRE-Zone_G3300.v7.5.inv.iso.incomp.visc"
#
curDir=${PWD##*/}
#
runfileTemplate=$3
PBStemplate=$4
run_out="$PWD/$5"
pbs_out=$6
#
cat $runfileTemplate | sed 's|NOD|'"$nod"'|' | sed 's|ELM|'"$elm"'|' | sed 's|BND|'"$bnd"'|' | sed 's|REGSTACK|'"$reg"'|' | sed 's|DSP|'"$dsp"'|' | sed 's|OUTPUT|'"$out"'|' > $run_out
#
#pbsName=$subj$(echo "$invDir" | tr '[:lower:]' '[:upper:]')
cat $PBStemplate | sed 's|DIR|'"$subj"'|' | sed 's|RUNFILE|'"$run_out"'|' > $pbs_out
#
qsub $pbs_out
