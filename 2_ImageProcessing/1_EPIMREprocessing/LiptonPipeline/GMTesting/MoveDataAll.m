mkdir('GraceData')
dirlist = dir('2023-*')

for ii=1:length(dirlist)
    mkdir(dirlist(ii).name)
end


for ii=1:length(dirlist)
    eval(sprintf('!cp %s/t2stack.mat GraceData/%s/t2stack.mat',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/t2stack.nii GraceData/%s/t2stack.nii',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/mre_for_inversion.mat GraceData/%s/mre_for_inversion.mat',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/t2bet.nii GraceData/%s/t2bet.nii',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/coT1W_3D_TFE.nii GraceData/%s/coT1W_3D_TFE.nii',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/NLI_Outputs/Mu.mat GraceData/%s/Mu.mat',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/NLI_Outputs/DR.mat GraceData/%s/DR.mat',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp %s/NLI_Outputs/Mu_%s.png GraceData/%s/Mu_%s.png',dirlist(ii).name,dirlist(ii).name,dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp -r %s/Nifti_Data GraceData/%s/Nifti_Data',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp -r %s/Registrations GraceData/%s/Registrations',dirlist(ii).name,dirlist(ii).name))
    eval(sprintf('!cp -r %s/FreeSurferOutput/mri/ GraceData/%s/FreeSurferOutput',dirlist(ii).name,dirlist(ii).name))

end