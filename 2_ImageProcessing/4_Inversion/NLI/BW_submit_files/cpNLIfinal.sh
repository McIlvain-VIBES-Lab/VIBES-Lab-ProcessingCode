#!/bin/bash
#
dirs=(EPILEPSY_110_mask070
EPILEPSY_110_mask080
EPILEPSY_110_mask090
EPILEPSY_110_mask095
EPILEPSY_111_mask070
EPILEPSY_111_mask080
EPILEPSY_111_mask090
EPILEPSY_111_mask095
EPILEPSY_112_mask070
EPILEPSY_112_mask080
EPILEPSY_112_mask090
EPILEPSY_112_mask095
EPILEPSY_113_mask070
EPILEPSY_113_mask080
EPILEPSY_113_mask090
EPILEPSY_113_mask095
EPILEPSY_114_mask070
EPILEPSY_114_mask080
EPILEPSY_114_mask090
EPILEPSY_114_mask095
EPILEPSY_116_mask070
EPILEPSY_116_mask080
EPILEPSY_116_mask090
EPILEPSY_116_mask095)
#
N=${#dirs[@]}
#
file='MRE-Zone_G3300.v7.5.inv.iso.incomp.visc..0100.ReconProps.mat'
file_paste=$file
dest_dir='/Volumes/Mediasonic/Curtis/160323_Epilepsy_completed/'
#
#inv='inv/'
inv='inv_SP_1e_09/'
for (( i=0; i<=$(( $N-1 )); i++ ))
do
  dir=${dirs[$i]}'/hex/'${dirs[$i]}'_voxelmesh/'$inv
  echo $dir
  mkdir -p $dest_dir$dir
  cp $dir$file $dest_dir$dir$file_paste
done
