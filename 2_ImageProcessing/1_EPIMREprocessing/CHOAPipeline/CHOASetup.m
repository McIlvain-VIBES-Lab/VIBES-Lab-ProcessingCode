function [info] =CHOASetup(subjpath)
disp('Setting Up Files')
disp('Do Not Change Your Location')
code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/CHOAPipeline';
addpath(code_path);

cd (subjpath); 
expected_root = '/Volumes/McIlvainDrive2/CHOA_Data/SUBJECT_DATA';
current_path = pwd;

[~, subj_id, ~] = fileparts(current_path);


if startsWith(current_path, expected_root)
    disp('The current path is correct.');
else
    error('This script must be run from within a subject folder inside: %s', expected_root);
end

disp('Changing Folder to study')
dir_list = dir; 
for i = 1:length(dir_list)
    if dir_list(i).isdir && contains(dir_list(i).name, 'M')
        movefile(dir_list(i).name, 'study');
        break; 
    end
end

disp('Pulling T1, MRE Phase and Mag out for study folder')
cd study/
cd (dir('Mri_*').name)

dir_list = dir;
for i = 1:length(dir_list)
    if strcmp(dir_list(i).name, '.') || strcmp(dir_list(i).name, '..')
        continue;
    end

    if dir_list(i).isdir && startsWith(dir_list(i).name, '3D_T1_MPRAGE_SAG')
        movefile(dir_list(i).name, ['../' filesep '3D_T1_MPRAGE_SAG']);
    end

    if dir_list(i).isdir && (contains(dir_list(i).name, 'XA60_Mag'))  
        movefile(dir_list(i).name, ['../' filesep 'MRE_3D_AX_ON_AXIS_Mag']);
    end
    if dir_list(i).isdir && startsWith(dir_list(i).name, 'XA60_P_P')
        movefile(dir_list(i).name, ['../' filesep 'MRE_3D_AX_ON_AXIS_P_P']);
    end
end
cd ..

disp('Changing Subject Folder Name')
dir1 = dir('MRE_3D_AX_ON_AXIS_Mag/*.dcm');
info = dicominfo(fullfile('MRE_3D_AX_ON_AXIS_Mag', dir1(1).name));
% new_folder = fullfile(expected_root, info.PatientID);
% movefile(current_path, new_folder);
% disp(new_folder)
% cd(new_folder)
% rmdir(current_path, 's');



