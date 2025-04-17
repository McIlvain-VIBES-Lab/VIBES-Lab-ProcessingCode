function LiptonCleanData(mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR)

%% Section 3 of the MRE Processing Code
% add path
addpath("/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/4_Inversion/NLI/NLIPrep/MREpreprocess/ToastPostProcess");

% Copy To NLI Folder
[~, SubjectName] = system('basename "$PWD"');
SubjectName = strtrim(SubjectName); 

load("mre_for_inversion.mat")
save(sprintf('%s.mat',SubjectName),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
% the subjectname.mat is saved (this is the mre for inversion mat)

UIUC_data_convert_mcilvain(SubjectName)
cd(sprintf('%s',SubjectName))
MRE_preprocess_v9_liang('default',SubjectName)
cd ..
eval(sprintf('!mv %s.mat /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))
eval(sprintf('!mv %s /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))

% full_path
full_path = "/Volumes/McIlvainDrive2/Send_to_NLI/2023-U7778-0290-MM-01";
% Ready for Testing (in liu of lines 34 and 35):
system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName)); 
% system('scp -r "/Volumes/McIlvainDrive2/Send_to_NLI/2023-U7778-0290-MM-01" gm3128@insomnia.rcs.columbia.edu:"/insomnia001/depts/mcilvain/users/mcilvain/2023-U7778-0290-MM-01"');

pause(30)
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))


end