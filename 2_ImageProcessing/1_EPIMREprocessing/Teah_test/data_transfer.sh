#!/bin/bash

subjects=(
'2023-U7487-0193-KC-01'
'2023-U7487-0377-KH-01'
'2023-U7487-0411-NP-01'
'2023-U7487-0424-TM-01'
'2023-U7487-0443-EM-01'
'2023-U7487-1040-LT-01'
'2023-U7487-1048-SL-01'
'2023-U7487-1052-CG-01'
'2023-U7487-1110-JA-01'
'2023-U7487-1136-IT-01'
'2023-U7487-1156-EM-01'
'2023-U7487-0989-MS-01'

)


for subject in "${subjects[@]}"; do
  cp -r "$subject/FreeSurferOutput" "/Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/FreeSurferOutput/$subject"
  #cp $subject/coT1W_3D_TFE.nii  /Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/t1_dicom/$subject.nii
  cp -r $subject/Nifti_Data /Volumes/McIlvainDrive2/Lipton_Soccer_Study/TeahTransfer/Nifti_Data/$subject
done

