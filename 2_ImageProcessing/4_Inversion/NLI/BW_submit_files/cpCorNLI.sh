#!/bin/bash
#
TemplateFiles='TemplateFiles/MRE-Zone_G3300.v7.5.inv.iso.incomp.visc.*'
#
dirs=(topup_mag_phs_00pos
topup_real_imag_00pos
topup_mag_phs_20pos
topup_real_imag_20pos
topup_mag_phs_20neg
topup_real_imag_20neg)
#
N=${#dirs[@]}
#
inv='inv/'
#inv='inv_SP_1e_09/'
#
for (( i=0; i<=$(( $N-1 )); i++ ))
do
  dir=${dirs[$i]}'/hex/'${dirs[$i]}'_voxelmesh/'$inv
  echo $dir
  cp $TemplateFiles $dir
done

