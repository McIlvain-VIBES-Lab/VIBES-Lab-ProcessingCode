#!/bin/bash

subjects=(
 '2023-U7778-00188-JG-01'
)


for subject in "${subjects[@]}"; do
  cp -r "$subject/FreeSurferOutput" "/Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/FreeSurferOutput/$subject"
  #cp $subject/coT1W_3D_TFE.nii  /Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/t1_dicom/$subject.nii
  cp -r $subject/Nifti_Data /Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/Nifti_Data/$subject
done

