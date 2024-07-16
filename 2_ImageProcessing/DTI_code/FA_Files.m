path='/Volumes/MRELAB2/EPILEPSY/';
subjects=dir([path,'1*']);
for i =1: length(subjects)
    copyfile([path, subjects(i).name, '/DTI_recon/dti_FA.nii'],...
        [path, 'FA_Files/',subjects(i).name,'_dti_FA.nii']);
end
