#!/bin/bash
#
# --------------------------
#
TemplateDir='TemplateFiles/'
pbs_template='MRE-Submit_template.pbs'
#
# Choose parameter set below
#
#runfileTemp='runfile-MRE-Zone.VE.50Hz.2p0mm.template'
#runfileOut='runfile-MRE-Zone.VE.50Hz.2p0mm.dat'
#pbsOut='MRE-Submit-v7p5_VE_50Hz_2p0mm.pbs'
#inv='inv'
#
#runfileTemp='runfile-MRE-Zone.VE.50Hz.2p5mm.template'
#runfileOut='runfile-MRE-Zone.VE.50Hz.2p5mm.dat'
#pbsOut='MRE-Submit-v7p5_VE_50Hz_2p5mm.pbs'
#inv='inv'
#
runfileTemp='runfile-MRE-Zone.VE.50Hz.3p0mm.template'
runfileOut='runfile-MRE-Zone.VE.50Hz.3p0mm.dat'
pbsOut='MRE-Submit-v7p5_VE_50Hz_3p0mm.pbs'
inv='inv'
#
# --------------------------
#
dirs=(topup_mag_phs_00pos
topup_real_imag_00pos
topup_mag_phs_20pos
topup_real_imag_20pos
topup_mag_phs_20neg
topup_real_imag_20neg)
#
N=${#dirs[@]}
startDir=`pwd`
#
for (( i=0; i<=$(( $N-1 )); i++ ))
do
  Qname=${dirs[$i]}
  #Qname=${dirs[$i]}'_SP'
  dir=${dirs[$i]}'/hex/'${dirs[$i]}'_voxelmesh/'
  cp $TemplateDir'build_submit_NLI.sh' $TemplateDir$pbs_template $TemplateDir$runfileTemp $dir
  cd $dir
  chmod +x build_submit_NLI.sh
  mkdir $inv
  ./build_submit_NLI.sh $Qname $inv $runfileTemp $pbs_template $runfileOut $pbsOut
  cd $startDir
done
