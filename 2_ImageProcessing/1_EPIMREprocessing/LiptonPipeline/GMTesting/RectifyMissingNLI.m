dirlist = dir('2023*')
for ii=1:length(dirlist)
 cd(dirlist(ii).name)
    if ~exist('nli_outputs','dir')

% Copy To NLI Folder
[~, SubjectName] = system('basename "$PWD"');
SubjectName = strtrim(SubjectName); 

save(sprintf('%s.mat',SubjectName),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
UIUC_data_convert_mcilvain(SubjectName)
cd(sprintf('%s',SubjectName))
MRE_preprocess_v9_mcilvain('default',SubjectName)
cd ..
eval(sprintf('!mv %s.mat /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))
eval(sprintf('!mv %s /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))

% Ready for Testing (in liu of lines 34 and 35):
system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName)); 
pause(20)
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))


    end 
    cd ..
end