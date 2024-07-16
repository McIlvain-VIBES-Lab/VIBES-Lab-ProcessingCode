%% Quicker Way to Do it, linear transform, fine for figures, not so good for pulling out values

!$FSLDIR/bin/flirt -ref MPRAGE_brain.nii -in t2bet.nii -out MRE_to_brain_flirt.nii -omat MRE_to_brain_flirt.mat
!$FSLDIR/bin/flirt -in Mu.nii -ref MPRAGE_brain.nii -out Mu2brain.nii -init MRE_to_brain_flirt.mat -applyxfm;
!$FSLDIR/bin/flirt -in DR.nii -ref MPRAGE_brain.nii -out DR2brain.nii -init MRE_to_brain_flirt.mat -applyxfm;


!$FSLDIR/bin/flirt -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -in MRE_to_brain_flirt.nii -out t2MNI.nii -omat t2MNI.mat;
!gunzip -f *.nii.gz
!$FSLDIR/bin/flirt -in Mu2brain.nii -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -out Mu2MNI.nii -init t2MNI.mat -applyxfm;
!$FSLDIR/bin/flirt -in DR2brain.nii -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -out DR2MNI.nii -init t2MNI.mat -applyxfm;

%% More robust way -- nonlinear reg, takes a little longer
!$FSLDIR/bin/bet2 MPRAGE.nii MPRAGE_brain.nii -m -v -f 0.25 -w 1.3
!$FSLDIR/bin/flirt -ref MPRAGE_brain.nii -in t2bet.nii -out MRE_to_brain_flirt.nii -omat MRE_to_brain_flirt.mat
!$FSLDIR/bin/flirt -in Mu.nii -ref MPRAGE_brain.nii -out Mu2brain.nii -init MRE_to_brain_flirt.mat -applyxfm;
!$FSLDIR/bin/flirt -in DR.nii -ref MPRAGE_brain.nii -out DR2brain.nii -init MRE_to_brain_flirt.mat -applyxfm;


!$FSLDIR/bin/flirt -in MRE_to_brain_flirt -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain -out MNI_to_MPR_flirt -omat MNI_to_MPR_flirt.mat -dof 12
!gunzip -f *.nii.gz
!$FSLDIR/bin/fnirt --ref=usr/local/fsl/data/standard/MNI152_T1_2mm_brain  --in=MPRAGE_brain.nii --aff=MNI_to_MPR_flirt.mat --cout=MNI_to_MPR_fnirt_cout.nii --iout=MNI_to_MPR_fnirt.nii
!$FSLDIR/bin/applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz --in=Mu2brain.nii --out=Mu2MNI_fnirt.nii --warp=brain_to_MNI_fnirt_cout.nii
!$FSLDIR/bin/applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz --in=DR2brain.nii --out=DR2MNI_fnirt.nii --warp=brain_to_MNI_fnirt_cout.nii
