function [info,SubjectName] = VarianSetup(subjpath,code_path)
disp('Setting Up Files')
disp('Do Not Change Your Location')
%code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
%addpath(code_path);

cd (subjpath); 
% TS 9/27/25
%Creates an Is Running File
mkdir('log')
%mkdir(fullfile(pwd, 'log')) %- Helen

expected_root = '/Volumes/McIlvainDrive2/Varian_Study/SUBJECT_DATA';
current_path = pwd;

[~, subj_id, ~] = fileparts(current_path);

logFilePath = fullfile('log/PreNLI_processing.txt');
fid = fopen(logFilePath, 'a');
if fid ~= -1
    fprintf(fid, 'Pre-NLI processing is Running for subject: %s\n', subj_id);
    fclose(fid);
end

if startsWith(current_path, expected_root)
    disp('The current path is correct.');
else
    error('This script must be run from within a subject folder inside: %s', expected_root);
end

disp('Changing Folder to study')
dir_list = dir; 
for i = 1:length(dir_list)
    if dir_list(i).isdir && startsWith(dir_list(i).name, 'A')
        movefile(dir_list(i).name, 'study');
        break; 
    end
end

disp('Pulling T1, Ax_Brain and QSM out for study folder')
cd study/
dir_list = dir;
for i = 1:length(dir_list)
    if strcmp(dir_list(i).name, '.') || strcmp(dir_list(i).name, '..')
        continue;
    end

for i = 1:length(dir_list)
    if strcmp(dir_list(i).name, '.') || strcmp(dir_list(i).name, '..')
        continue;
    end

    if dir_list(i).isdir && startsWith(dir_list(i).name, 'Ax_BRAIN_MRE')
        try
            movefile(dir_list(i).name, ['..' filesep 'Ax_BRAIN_MRE']);
        catch ME
            fprintf('Could not move folder %s: %s\n', dir_list(i).name, ME.message);
            continue;
        end
    end

    if dir_list(i).isdir && (startsWith(dir_list(i).name, 'SAGITTAL_T1_VOLUME') ) %|| startsWith(dir_list(i).name, 'ORIG_102_T1W_3D_TFE'))
        try
            movefile(dir_list(i).name, ['..' filesep 'SAGITTAL_T1_VOLUME']);
        catch ME
            fprintf('Could not move folder %s: %s\n', dir_list(i).name, ME.message);
            continue;
        end
    end

end

end
cd ..

mkdir('Archive');
mkdir('Test_Dev');

disp('Changing Subject Folder Name')
dir1 = dir('Ax_BRAIN_MRE/*.dcm');
info = dicominfo(fullfile('Ax_BRAIN_MRE', dir1(1).name));
SubjectName = info.PatientID;



