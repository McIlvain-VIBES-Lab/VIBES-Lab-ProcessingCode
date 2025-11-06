%% Lipton Wrapper Code
% October 10th 2025
% Grace McIlvain Teah Serani Hailin Liang
% Helen's New Modification Explain:
% Prepare multi-masks 
% loop through all masks (whether multi predictions / multi human masks) 
% Send to NLI simultaneously
% 1. Requires edits on the HelenLiptonCleanData (run only after all masks
% send to NLI)
% 2. Requires edits on subject folder name to prevent being overwritten
% after starting NLI solver (insomnia)

% Pls dont delete this code as we may implement this code for Helen's
% future project


%% Start up Code (No Change) 
clear all; close all;
delete(gcp('nocreate'));
disp('Code Starting Please Wait')

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = ['/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/'];
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general


%% SETUP YOUR SUBJECT NAME
% THIS STEP IS CRUCIAL , plssss run 
pwd;
SubjectName = pwd;

%% Modify Teah's Lipton Soccer Seteup Function
disp('Setting Up Files')
disp('Do Not Change Your Location')
%code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
%addpath(code_path);
subjpath = pwd; 
cd (subjpath); 
% TS 9/27/25
%Creates an Is Running File
mkdir('log')
% mkdir(fullfile(pwd, 'log')) - Helen

expected_root = '/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA';
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

    if dir_list(i).isdir && (startsWith(dir_list(i).name, '102_T1W_3D_TFE') ) % || startsWith(dir_list(i).name, 'ORIG_102_T1W_3D_TFE'))
        try
            movefile(dir_list(i).name, ['..' filesep 'T1W_3D_TFE']);
        catch ME
            fprintf('Could not move folder %s: %s\n', dir_list(i).name, ME.message);
            continue;
        end
    end

    if dir_list(i).isdir && startsWith(dir_list(i).name, '103_3D_Ax_QSM')
        try
            movefile(dir_list(i).name, ['..' filesep 'QSM']);
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
new_folder = fullfile(expected_root, SubjectName);
movefile(current_path, new_folder);
disp(new_folder)
cd(new_folder)
rmdir(current_path, 's');
% ignore error msg for now, move on!

%% run(fullfile(startup_path, 'startup_matlab_general.m'));

% Now you can call the LiptonSoccerSetup function from anywhere
% function added by TS March 26th 2025
% This will rename the subject dir and pull&name T1W and Ax_Brain
% subjpath = pwd; 
% [info,SubjectName] = LiptonSoccerSetup(subjpath,code_path);

% Section 1
dx = (info.ReconstructionDiameter)/double(info.Height);
dy = (info.ReconstructionDiameter)/double(info.Width);
dz = info.SliceThickness; 
freq =double(50); % TS: We need to check where the frequency is in the header

cd('Ax_BRAIN_MRE')
[phsimg,t2stack] = GE_MRE_ProcessingPart1(dx,dy,dz,freq); 

% Anatomical Scan DCM to NII
%cd('T1W_3D_TFE')
T1_dcm2nii 
%cd ..


%% Helen's Modification 
% We will like to run the section 2 "with our experts' manual masking 
% need to redefine the variable "mask" which will go to the second
% Processing
% below set contains all the masks
% depends on your needs, the maskx could be the pure t2bet mask or the
% multiple predictions/multi experts' masks

% below are the maskx (altho it named 'mask')
t2bet = load('t2mask_bet.mat');
all_masks = {'G012_AH_mask.mat', 'G012_GM_mask.mat','G012_ZS_mask.mat','G012_TS_mask.mat'};
for i = 1:length(all_masks)
    mask_file = all_masks{i};

    % Extract initials
    tokens = split(mask_file, {'_', '.'});
    if numel(tokens) >= 2
        mask_label = tokens{2};
    else
        error('Unexpected mask file format: %s', mask_file);
    end

    % Load mask
    % using the key "maskx"
    loaded = load(mask_file);
    if isfield(loaded, 'maskx')
        manual_mask = loaded.maskx;

    else
        error('The file %s does not contain a variable named "maskx"', mask_file);
    end

    fprintf('\n=== Processing with manual mask: %s ===\n', mask_label);
   

    % Section 2: Processing
    [mreParams, mask, Zmotion, Ymotion, Xmotion, t2stack, OSS_SNR] = ...
       Helen_GE_MRE_ProcessingPart2(phsimg, t2stack, t2bet, manual_mask, dx, dy, dz, freq, mask_label);

    % Section 3A: Send to NLI (mask version)
    HelenLiptonPrepareForNLI(mreParams, mask, Zmotion, Ymotion, Xmotion, t2stack, OSS_SNR, mask_label);
end

% Section 3B: Final cleanup (once all sent)
% HelenLiptonCleanData();


%% Section 2
% !mkdir dcfiles
% [mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR] = GE_MRE_ProcessingPart2(phsimg,t2stack,mask,dx,dy,dz,freq);

%% Section 3
% clean up data
% send things to NLI??? Helen
% Helen adds a mask_label param in the HelenLiptionCleanData Below
% HelenLiptonCleanData(mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR)

%% Section 4
%disp("Running Freesurfer on Insomnia")
%system(['sh /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline/freesurfer_insomnia.sh']);

%% Section 5 
% TS 9/27/25
% % Delete the IsRunning File and creates a complete file 
% delete(fullfile('log', 'PreNLI_processing.txt'));
% 
% logFilePath = fullfile('log/PreNLI_complete.txt');
% fid = fopen(logFilePath, 'a');
% if fid ~= -1
%     fprintf(fid, 'Pre-NLI processing complete for subject: %s\n', SubjectName);
%     fclose(fid);
% end
% 
% logFilePath = fullfile('log/FS_processing.txt');
% fid = fopen(logFilePath, 'a');
% if fid ~= -1
%     fprintf(fid, 'FreeSurfer running for subject: %s\n', SubjectName);
%     fclose(fid);
% end
% 
% %TS Needs to move Section 5 
% %% Section 5
% %disp("Running Registerations on Insomnia")
% %system(['sh /Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline/reg_insomnia.sh'])