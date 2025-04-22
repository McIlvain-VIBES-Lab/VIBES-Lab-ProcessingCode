function LiptonCleanData(mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR)

%% Section 3 of the MRE Processing Code
% Clean Lipton Data 
% March 26th 2025
% Grace McIlvain
% Clean Up data

mkdir('File_Storage')
!mv Ax_Brain_MRE/ File_Storage
!mv study/ File_Storage
!mv T1W_3D_TFE/ File_Storage/
!mv dcfiles/ File_Storage/

!mv maskx.mat File_Storage/
!mv mre_mag.nii File_Storage/
!mv mre_phs.nii File_Storage/
!mv mre_mask.nii File_Storage/
!mv mre_output.nii File_Storage/

!mv mreimages_unwrap.mat File_Storage/
!mv t2mask_bet.mat File_Storage/
% 
% 
% Copy To NLI Folder
[~, SubjectName] = system('basename "$PWD"');

SubjectName = strtrim(SubjectName); 
addpath(SubjectName)

save(sprintf('%s.mat',SubjectName),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
UIUC_data_convert_mcilvain(SubjectName)
cd(sprintf('%s',SubjectName))
MRE_preprocess_v9_mcilvain('default',SubjectName)
cd ..
%eval(sprintf('!mv %s.mat /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))
%eval(sprintf('!mv %s /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))

% Ready for Testing (in liu of lines 34 and 35):
system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName)); 
pause(20)
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))



end