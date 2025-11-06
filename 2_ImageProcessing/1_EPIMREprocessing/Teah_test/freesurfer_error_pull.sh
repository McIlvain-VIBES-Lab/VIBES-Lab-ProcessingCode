#!/bin/bash

# Set local path and directories to check
local_path="/Volumes/McIlvainDrive2/Lipton_Soccer_Study/sandbox"
remote_path="/insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs"

# Get list of subject directories (2023-U*)
dirs=($(find "$local_path" -maxdepth 1 -type d -name '2023-U*' -exec basename {} \;))

# Print the directories for debugging
echo "Directories to check:"
for dir in "${dirs[@]}"; do
    echo "$dir"
done
echo "-----------------------------------------------------"
# Initialize an array to track subjects needing transfer
subject_list=()

# Loop through subject directories and check for missing FreeSurferOutput
for subj in "${dirs[@]}"; do
    echo "Checking local directory: $subj"

    if [ -d "$local_path/$subj/FreeSurferOutput" ]; then
        echo "FreeSurfer Output already on local machine"
    else
        subject_list+=("$subj")
    fi
done

# Print out the subjects missing FreeSurferOutput
echo "Subjects missing FreeSurfer Output:"
for s in "${subject_list[@]}"; do
    echo "$s"
done
echo "-----------------------------------------------------"
echo "Subjects with FreeSurfer Outputs on Insomnia:"
subject_list_str="${subject_list[@]}"

existing_subjects=$(ssh ts3641@insomnia.rcs.columbia.edu bash <<EOF
for subject in $subject_list_str; do
    if [ -f "$remote_path/\$subject/scripts/recon-all.done" ] && [ ! -f "$remote_path/\$subject/scripts/recon-all.error" ]; then
        echo \$subject
    fi
done
EOF
)

new_subjects=$(for exist in $existing_subjects; do
    if [ ! -d "$local_path/$exist/FreeSurferOutput" ]; then
        echo "$exist"
    fi
done)

echo $new_subjects

if [ -z "$new_subjects" ]; then
    echo "All subjects FreeSurfer outputs are on local computer. Nothing to do."
    exit 0
fi

# Prepare list of remote paths to rsync
rsync_sources=()
for subject in $new_subjects; do
    rsync_sources+=("ts3641@insomnia.rcs.columbia.edu:$remote_path/$subject")
done

temp_download_subject="$local_path/TEMP_TRANSFER_DATA"
mkdir -p "$temp_download_subject"

echo "Starting rsync of the following subjects:"
echo "$new_subjects"

rsync -avz --relative --progress "${rsync_sources[@]}" "$temp_download_subject"

for subject in $new_subjects; do
    src="$temp_download_subject/insomnia001/depts/mcilvain/users/mcilvain/FreeSurferOutputs/$subject"
    dest="$local_path/$subject/FreeSurferOutput"

    echo "Moving $subject -> $dest"
    mkdir -p "$local_path/$subject"

    if [ -d "$src" ]; then
        mv "$src" "$dest"
    else
        echo "ERROR: Source directory $src does not exist, skipping move."
    fi
done

