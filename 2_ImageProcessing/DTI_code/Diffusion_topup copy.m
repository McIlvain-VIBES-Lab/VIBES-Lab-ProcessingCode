function Diffusion_topup

%% Must be run from Diffusion folder that has a topup folder in it
%!$FSLDIR/bin/topup 

%cd Buckley_Tbi_Brain/Diffusion
mkdir topup
movefile Diff_AP/Diff_AP.nii Diff_AP.nii 
movefile Diff_PA/Diff_PA.nii Diff_PA.nii  
movefile acq_params_topup.txt topup/acq_params_topup.txt
!$FSLDIR/bin/fslroi Diff_AP.nii Diff_AP_B0.nii 0 1
!$FSLDIR/bin/fslroi Diff_PA.nii Diff_PA_B0.nii 0 1
!gunzip -f *.nii.gz
!$FSLDIR/bin/fslmerge -t topup/mergeAP_PA Diff_AP_B0.nii Diff_PA_B0.nii
cd topup
!gunzip -f *.nii.gz

%start with location
% --imain=(filename)
!$FSLDIR/bin/topup --imain=mergeAP_PA.nii --datain=acq_params_topup.txt --iout=topup_field --out=topup_results -v
cd ..
fprintf('Applying Topup')
!$FSLDIR/bin/applytopup --imain=Diff_AP.nii,Diff_PA.nii --datain=topup/acq_params_topup.txt --inindex=1,2  --topup=topup/topup_results --out=topup/Diff_topup -v
cd topup
!gunzip -f *.nii.gz
DTI_img = load_untouch_nii('Diff_topup.nii');
DTI = DTI_img.img; 
cd ..; save DTI.mat DTI

%% Reorganize Folders
movefile topup/Diff_topup.nii Diff_topup.nii 
mkdir DiffCreation
dirval = dir('Diff_AP/*.bval'); dirvec= dir('Diff_AP/*.bvec');
movefile(sprintf('Diff_AP/%s',dirval.name),'Diff.bval')
movefile(sprintf('Diff_AP/%s',dirvec.name),'Diff.bvec')
movefile Diff_AP.nii DiffCreation/Diff_AP.nii
movefile Diff_AP_B0.nii DiffCreation/Diff_AP_B0.nii
dirval = dir('Diff_PA/*.bval'); dirvec= dir('Diff_PA/*.bvec');
movefile(sprintf('Diff_PA/%s',dirval.name),'DiffCreation/Diff.bval')
movefile(sprintf('Diff_PA/%s',dirvec.name),'DiffCreation/Diff.bvec')
movefile Diff_PA.nii DiffCreation/Diff_PA.nii
movefile Diff_PA_B0.nii DiffCreation/Diff_PA_B0.nii

% --imain
% --datain
% --iout 