#!/bin/bash

subjects=(
"2023-U7778-0533-PD-01"
"2023-U7778-0735-DF-01"
"2023-U7778-0779-PG-01"

)


# )

for subject in "${subjects[@]}"; do
  
  #sh /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/Teah_test/reg_insomnia.sh 

  # mkdir /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject

  # cp /Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/$subject/NLI_Outputs/Mu.mat /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Mu.mat 
  # cp /Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/$subject/Registrations/aparc+aseg-reg.nii /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/aparc+aseg-reg.nii
  # cp /Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/$subject/Registrations/brain-reg.nii /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/brain-reg.nii
  # cp /Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/$subject/Registrations/Magnitude.nii.gz /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Magnitude.nii.gz


  # rsync -avz "/Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA/$subject/FreeSurferOutput/mri/ribbon.mgz" "ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/Registrations/$subject/ribbon.mgz" 

  rsync -avz "ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/Registrations/$subject/ribbon-reg.nii" "/Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii"
  
  # mkdir /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject
  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Left-Cerebral-White-Matter.nii --match 2
  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Left-Cerebral-Cortex.nii --match 3

  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Right-Cerebral-White-Matter.nii --match 41
  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Right-Cerebral-Cortex.nii --match 42

  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Cerebral-White-Matter.nii --match 41 --match 2
  mri_binarize --i /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/ribbon-reg.nii --o /Volumes/McIlvainDrive2/Lipton_Soccer_Study/GraceSegmentations/$subject/Cerebral-Cortex.nii --match 42 --match 3

done
