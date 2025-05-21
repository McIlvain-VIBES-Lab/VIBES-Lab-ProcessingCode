#!/bin/bash

# Get the current directory and subject ID
subj_dir=$(pwd)
subj_id=$(basename $(pwd))

# Define the file system directory on the remote server
fs_dir="ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs"

# Check if the 'coT1W_3D_TFE.nii' file exists
if [ -f "$subj_dir/coT1W_3D_TFE.nii" ]; then
    # Copy the file to the remote server
    scp "$subj_dir/coT1W_3D_TFE.nii" $fs_dir/$subj_id.nii
    
    # SSH into the server and submit the job
    ssh ts3641@insomnia.rcs.columbia.edu <<EOF
    cd /insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs
    sbatch lipton_fs.sh $subj_id.nii
EOF

# Check if the 'T1W_3D_TFE.nii' file exists
elif [ -f "$subj_dir/T1W_3D_TFE.nii" ]; then
    # Copy the file to the remote server
    scp "$subj_dir/T1W_3D_TFE.nii" $fs_dir/$subj_id.nii
    
    # SSH into the server and submit the job
    ssh ts3641@insomnia.rcs.columbia.edu <<EOF
    cd /insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs
    sbatch lipton_fs.sh $subj_id.nii
EOF

# If neither file is found, print an error message
else
    echo "Subject missing T1"
fi
