%% HCP Diff Code processing code

dirs = dir('Aniso*');
for ii = 1:length(dirs)-2
        cd(dirs(ii).name)
        cd Cljohnson_R01_Anisotropy/Diffusion/HCPdiff
        cd HCPdiff
        !rm *.nii
        !rm *.mat
        !/Applications/MRIcron/dcm2nii * -4fn
        !gunzip -f *.nii.gz
        dirnii = dir('*.nii'); movefile(sprintf('%s',dirnii.name),'DTI.nii')
        dirbvec = dir('*.bvec'); movefile(sprintf('%s',dirbvec.name),'Diff.bvec')
        dirbval = dir('*.bval'); movefile(sprintf('%s',dirbval.name),'Diff.bval')
        
        cd ../Diff_AP
        !rm *.nii
        !rm *.mat
        !/Applications/MRIcron/dcm2nii * -4fn
        !gunzip -f *.nii.gz
        dirnii = dir('*.nii'); movefile(sprintf('%s',dirnii.name),'Diff_AP.nii')
        
        cd ../Diff_PA
        !rm *.nii
        !rm *.mat
        !/Applications/MRIcron/dcm2nii * -4fn
        !gunzip -f *.nii.gz
        dirnii = dir('*.nii'); movefile(sprintf('%s',dirnii.name),'Diff_PA.nii')
        cd ..
        
        Diffusion_topup_HCPdiff
        dti_proc_topup_HCPdiff
        mkdir Diff_rotate
        !mv t2bet.nii Diff_rotate/t2bet.nii
        !mv dti_S0.nii Diff_rotate/dti_S0.nii
        !mv Diff_topup.nii Diff_rotate/Diff_topup.nii
        !mv Diff.bvec Diff_rotate/Diff.bvec
        !mv Diff.bval Diff_rotate/Diff.bval
        cd Diff_rotate
        dti_proc_rotatebvecs
        cd ..
        cd ../../../..
end