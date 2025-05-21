#!/bin/bash

# Your local path and directories to check
local_path="/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA/"

# Find directories with the pattern and extract basenames
dirs=($(find "$local_path" -maxdepth 1 -type d -name '2023-U*' -exec basename {} \;))

# Remote path on HPC
remote_path="/insomnia001/depts/mcilvain/users/mcilvain/Registrations"

echo "Directories to check:"
for dir in "${dirs[@]}"; do
    echo "$dir"
    
    # Construct remote directory path
    remote_dir="$remote_path/$dir"
    # SSH and check if the 'Registration' file exists in the remote directory
    check=$(ssh -v gm3128@insomnia.rcs.columbia.edu bash <<EOF
    if [ -d "$remote_dir" ]; then
        echo 1
    else
        echo 0
    fi
EOF
)
    check=1
    echo $check

    if [ "$check" == "1" ]; then
        echo "$dir: Registrated without errors. Transferring Registations"
        echo "$local_path/$dir/Registrations"
        
        scp -r gm3128@insomnia.rcs.columbia.edu:"$remote_dir" "$local_path/$dir/Registrations"
    elif [ "$check" == "0" ]; then
        echo "$dir: Registrations does not exist on Insomnia"
    else
        echo "$dir: Unknown issue (check='$check')"
    fi
done