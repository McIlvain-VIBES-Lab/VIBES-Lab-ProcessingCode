#!/bin/bash

pwd
subj_dir=$(pwd)
subj_id=$(basename "$subj_dir")

fs_dir="/insomnia001/depts/mcilvain/users/mcilvain/Registrations/$subj_id"

ssh ts3641@insomnia.rcs.columbia.edu "mkdir -p $fs_dir"

if [ -f "$subj_dir/coT1W_3D_TFE.nii" ] && [ -f "$subj_dir/Nifti_Data/Magnitude.nii.gz" ] && [ -f "$subj_dir/Nifti_Data/Stiffness.nii.gz" ]; then
    scp "$subj_dir/coT1W_3D_TFE.nii" "$subj_dir/Nifti_Data/Magnitude.nii.gz" "$subj_dir/Nifti_Data/Stiffness.nii.gz" ts3641@insomnia.rcs.columbia.edu:$fs_dir
    ssh ts3641@insomnia.rcs.columbia.edu "cd $fs_dir && sbatch /insomnia001/depts/mcilvain/users/mcilvain/Registrations/lipton_reg.sh $subj_id"

else
    echo "Missing T1, t2stack.nii, or Stiffness.nii"
fi
