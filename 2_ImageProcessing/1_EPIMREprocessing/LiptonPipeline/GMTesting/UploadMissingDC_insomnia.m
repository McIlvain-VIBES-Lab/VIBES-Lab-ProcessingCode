dirlist = dir('2023-*')

for ii=1:length(dirlist)
    cd(dirist(ii).name)

    [~, SubjectName] = system('basename "$PWD"');
cd('FileStorage')
mkd
SubjectName = strtrim(SubjectName); 
disp(['SubjectName = ', SubjectName])

addpath(fullfile('/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA', SubjectName))


eval(sprintf('!mv %s.mat %s/',SubjectName,SubjectName))
cd ..


% Ready for Testing (in liu of lines 34 and 35):
system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName)); 
pause(20)
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
%system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))

