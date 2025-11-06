#!/bin/bash

pwd
subj_dir=$(pwd)
mkdir Registrations
subj_id=$(basename "$subj_dir")
echo $subj_id
fs_dir="/insomnia001/depts/mcilvain/users/mcilvain/Registrations/$subj_id"

if [ -f "$subj_dir/coT1W_3D_TFE.nii" ] && \
   [ -f "$subj_dir/Nifti_Data/Magnitude.nii.gz" ] && \
   [ -f "$subj_dir/Nifti_Data/Stiffness.nii.gz" ] && \
   [ -f "$subj_dir/FreeSurferOutput/mri/brain.mgz" ] && \
   [ -f "$subj_dir/FreeSurferOutput/mri/nu.mgz" ] && \
   [ -f "$subj_dir/FreeSurferOutput/mri/aparc+aseg.mgz" ]; then

    cp "$subj_dir/coT1W_3D_TFE.nii" "$subj_dir/Registrations/coT1W_3D_TFE.nii"
    cp "$subj_dir/Nifti_Data/Magnitude.nii.gz" "$subj_dir/Registrations/Magnitude.nii.gz"
    cp "$subj_dir/Nifti_Data/Stiffness.nii.gz" "$subj_dir/Registrations/Stiffness.nii.gz"
    cp "$subj_dir/FreeSurferOutput/mri/brain.mgz" "$subj_dir/Registrations/brain.mgz"
    cp "$subj_dir/FreeSurferOutput/mri/nu.mgz" "$subj_dir/Registrations/nu.mgz"
    cp "$subj_dir/FreeSurferOutput/mri/aparc+aseg.mgz" "$subj_dir/Registrations/aparc+aseg.mgz"

    rsync -avz "$subj_dir/Registrations/" "ts3641@insomnia.rcs.columbia.edu:$fs_dir/"
    ssh ts3641@insomnia.rcs.columbia.edu "cd $fs_dir && sbatch /insomnia001/depts/mcilvain/users/mcilvain/Registrations/lipton_reg.sh $subj_id"

else
    echo "Missing T1, t2stack.nii, Stiffness.nii or FreeSurfer Output"
fi
