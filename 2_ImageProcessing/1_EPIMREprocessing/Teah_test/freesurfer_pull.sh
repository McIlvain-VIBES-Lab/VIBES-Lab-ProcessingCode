#!/bin/bash

# Your local path and directories to check
local_path="/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA"
dirs=($(find "$local_path" -maxdepth 1 -type d -name '2023-U*' -exec basename {} \;))

# Remote path on HPC
remote_path="/insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs"

# Print the directories (for debugging purposes)
echo "Directories to check:"
for dir in "${dirs[@]}"; do
    echo "$dir"
done

for dir in "${dirs[@]}"; do
    echo "Checking directory: $dir"

    if [ -d "$local_path/$dir/FreeSurferOutput" ]; then
        echo "FreeSurfer Output already on local machine"
        continue
    fi

    remote_dir="$remote_path/$dir"

    check=$(ssh -T ts3641@insomnia.rcs.columbia.edu bash <<EOF
        if [ -d "$remote_dir" ]; then
            if [ -f "$remote_dir/scripts/recon-all.done" ] && [ ! -f "$remote_dir/scripts/recon-all.error" ]; then
                echo 1
            elif [ -f "$remote_dir/scripts/recon-all.error" ]; then
                echo 2
            elif [ -f "$remote_dir/scripts/IsRunning.lh+rh" ] && find "$remote_dir/scripts/IsRunning.lh+rh" -mmin +720 -print -quit | grep -q .; then
                echo 3
                rm -f "$remote_dir/scripts/IsRunning.lh+rh"
                sbatch "$remote_dir/lipton_rerun_fs.sh" "$dir"
            else
                echo 0
            fi
        else
            echo -1
        fi
EOF
    )

    if [ "$check" == "1" ]; then
        echo "$dir: Completed without errors. Transferring Outputs"
        scp -r ts3641@insomnia.rcs.columbia.edu:"$remote_dir" "$local_path/$dir/FreeSurferOutput"
    elif [ "$check" == "2" ]; then
        echo "$dir: Completed with errors. Transferring Error File"
        scp -r ts3641@insomnia.rcs.columbia.edu:"$remote_dir/scripts/recon-all.error" "$local_path/$dir/FreeSurfer.error"
    elif [ "$check" == "3" ]; then
        echo "$dir: Rerunning Recon-all"
    elif [ "$check" == "0" ]; then
        echo "$dir: Freesurfer still is Running"
    elif [ "$check" == "-1" ]; then
        echo "$dir: Directory does not exist on remote"
    else
        echo "$dir: Unknown issue (check='$check')"
    fi
done
