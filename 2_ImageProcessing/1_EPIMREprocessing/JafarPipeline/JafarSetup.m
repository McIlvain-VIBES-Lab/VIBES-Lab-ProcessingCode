function JafarSetup(subjpath)
disp('Setting Up Files')
disp('Do Not Change Your Location')
code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/JafarPipeline';
addpath(code_path);

cd (subjpath); 
expected_root = '/Volumes/McIlvainDrive2/TestandDev/Jafar/Dead_Brains';
current_path = pwd;

[~, subj_id, ~] = fileparts(current_path);


if startsWith(current_path, expected_root)
    disp('The current path is correct.');
else
    error('This script must be run from within a subject folder inside: %s', expected_root);
end

% disp('Changing Folder to study')
% dir_list = dir; 
% for i = 1:length(dir_list)
%     if dir_list(i).isdir && startsWith(dir_list(i).name, 'A')
%         movefile(dir_list(i).name, 'study');
%         break; 
%     end
% end
% 
% disp('Pulling T1, Ax_Brain and QSM out for study folder')
% cd study/
% dir_list = dir;
% for i = 1:length(dir_list)
%     if strcmp(dir_list(i).name, '.') || strcmp(dir_list(i).name, '..')
%         continue;
%     end
% 
%     if dir_list(i).isdir && startsWith(dir_list(i).name, 'Ax_BRAIN_MRE')
%         %movefile(dir_list(i).name, ['..' filesep 'Ax_BRAIN_MRE']);
%         eval(sprintf('!cp -r %s ../%s',dir_list(i).name,dir_list(i).name))
%     end
% 
%     % if dir_list(i).isdir && (startsWith(dir_list(i).name, '102_T1W_3D_TFE') || startsWith(dir_list(i).name, 'ORIG_102_T1W_3D_TFE'))
%     %     movefile(dir_list(i).name, ['..' filesep 'T1W_3D_TFE']);
%     % end
%     % if dir_list(i).isdir && startsWith(dir_list(i).name, 'QSM')
%     %     movefile(dir_list(i).name, ['..' filesep 'QSM']);
%     % end
% end
%cd ..

dirlist2= dir('Ax_BRAIN*')
for jj=1:length(dirlist2)
    cd(dirlist2(jj).name)
    mkdir('dicoms')
    addpath('dicoms')
    !mv *.dcm dicoms/
    
    cd ..
end


%disp('Changing Subject Folder Name')
%
%
% new_folder = fullfile(expected_root, info.PatientID);
% movefile(current_path, new_folder);
% disp(new_folder)
% cd(new_folder)
% rmdir(current_path, 's');


