#!/bin/bash

#LIPTON NLI PULL 
#This will check insomnia for complete NLI scripts and pull them back to local machine
#This script will ssh into insomnia once and pull NLI Outputs that are not on local machine
#Designed to run in checkfornl.m 

remote_path="/insomnia001/depts/mcilvain/users/mcilvain"


 if [ "$(pwd)" == "/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA" ]; then
 	local_path="/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA"
 fi
 if [ "$(pwd)" == "/Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA" ]; then
 	local_path="/Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA"
 fi


temp_download_subject="$local_path/TEMP_NLI_TRANSFER_DATA"

mkdir "$temp_download_subject"

# Get list of local subjects matching pattern
subjects=($(find "$local_path" -maxdepth 1 -type d -name '2023-U*' -exec basename {} \;))

# Convert array to space-separated string
subject_list="${subjects[@]}"
echo $subject
echo "Checking remote directories..."

existing_subjects=$(ssh ts3641@insomnia.rcs.columbia.edu bash -s -- "${subjects[@]}" <<'EOF'
for subject in "$@"; do
    [ -z "$subject" ] && continue
    subject_dir="/insomnia001/depts/mcilvain/users/mcilvain/${subject}/hex/${subject}_voxelmesh/inv/"
    if find "$subject_dir" -type f -name '*0100.prop*' 2>/dev/null | grep -q .; then
        echo "$subject"
    fi
done
EOF
)


# Filter to get only the directories that don't already exist locally
new_subjects=$(for exist in $existing_subjects; do 
    if [ ! -d "$local_path/$exist/NLI_Outputs" ]; then
        echo "$exist"
    fi
done)

echo $new_subjects
# # If no new subjects, exit early
# if [ -z "$new_subjects" ]; then
#     echo "All subjects registration files are on local computer. Nothing to do."
#     #exit 0
# fi

# Prepare list of remote paths to rsync
rsync_sources=()
for subject in $new_subjects; do
    echo "hi"
    rsync_sources+=("ts3641@insomnia.rcs.columbia.edu:$remote_path/$subject")
done

echo "Starting rsync of the following subjects:"
echo "$new_subjects"

# Use array expansion to rsync
rsync -avz --relative "${rsync_sources[@]}" "$temp_download_subject"

# Move them to the correct locations under each subject subject
for subject in $new_subjects; do
    src="$temp_download_subject/insomnia001/depts/mcilvain/users/mcilvain/$subject"
    dest="$local_path/$subject/NLI_Outputs"
    echo "Moving $subject -> $dest"
    mkdir -p "$local_path/$subject"  
    mv "$src" "$dest"
    echo "NLI output has returned for subject: $subject" >> "$local_path/$subject/log/needs_postNLI.txt"
done

# Clean up
rm -rf "$temp_download_subject"
echo "Transfer complete. Temp data removed."