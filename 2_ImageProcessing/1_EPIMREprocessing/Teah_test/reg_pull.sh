#!/bin/bash

# Local base path
local_path="/Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA"

# Remote base path
remote_path="/insomnia001/depts/mcilvain/users/mcilvain/Registrations"

# Temp download subjectectory
temp_download_subject="$local_path/TEMP_TRANSFER_DATA"

# Clean temp subjectectory

mksubject "$temp_download_subject"
 
# Get list of local subjects matching pattern
subjects=($(find "$local_path" -maxdepth 1 -type d -name '2023-U*' -exec basename {} \;))

# Convert array to space-separated string
subject_list="${subjects[@]}"

echo "Checking remote directories..."
# SSH once to check which subjects exist remotely
existing_subjects=$(ssh ts3641@insomnia.rcs.columbia.edu bash <<EOF
for subject in $subject_list; do
    if [ -d "$remote_path/\$subject" ]; then
        echo \$subject
    fi
done
EOF
)

#Filter to get only the directories that don't already exist locally
new_subjects=$(for exist in $existing_subjects; do 
    if [ ! -d "$local_path/$exist/Registrations" ]; then
        echo "$exist"
    fi
done)

# If no new subjects, exit early
if [ -z "$new_subjects" ]; then
    echo "All subjects registration files are on local computer. Nothing to do."
    exit 0
fi

# Prepare list of remote paths to rsync
rsync_sources=()
for subject in $new_subjects; do
    rsync_sources+=("ts3641@insomnia.rcs.columbia.edu:$remote_path/$subject")
done

echo "Starting rsync of the following subjects:"
echo "$new_subjects"

# Use array expansion to rsync
rsync -avz --relative "${rsync_sources[@]}" "$temp_download_subject"

# Move them to the correct locations under each subject subject
for subject in $new_subjects; do
    src="$temp_download_subject/insomnia001/depts/mcilvain/users/mcilvain/Registrations/$subject"
    dest="$local_path/$subject/Registrations"

    echo "Moving $subject -> $dest"
    mksubject -p "$local_path/$subject"
    mv "$src" "$dest"

    echo "converting .lta files to .mat files"
    lta_convert \
          --inlta $local_path/$subject/Registrations/aff_t1_to_t2.lta \
          --outfsl $local_path/$subject/Registrations/aff_t1_to_t2.mat \
          --src $local_path/$subject/Registrations/coT1W_3D_TFE.nii \
          --trg $local_path/$subject/Registrations/Magnitude.nii.gz

        lta_convert \
           --inlta $local_path/$subject/Registrations/aff_t2_to_t1.lta \
           --outfsl $local_path/$subject/Registrations/aff_t2_to_t1.mat \
           --src $local_path/$subject/Registrations/Magnitude.nii.gz \
           --trg $local_path/$subject/Registrations/coT1W_3D_TFE.nii

done

# Clean up
rm -rf "$temp_download_subject"
echo "Transfer complete. Temp data removed."