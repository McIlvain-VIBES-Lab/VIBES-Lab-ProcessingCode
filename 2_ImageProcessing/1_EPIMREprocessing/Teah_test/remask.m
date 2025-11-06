%% Lipton Remasking Code
% Oct 4, 2025
% Allen Hong Grace McIlvain Teah Serani

clear all; close all;

disp('Code Starting Please Wait')

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general


% TS added 11/4/25 for log files and sending pulling back from insomnia
delete('log/PostNLI_complete.txt')
delete('log/PreNLI_complete.txt')

[~, SubjectName] = system('basename "$PWD"');

logFilePath = fullfile('log/Remasking.txt');
fid = fopen(logFilePath, 'a');
if fid ~= -1
    fprintf(fid, 'Remasking for  subject: %s\n', SubjectName);
    fclose(fid);
end


SubjectName = strtrim(SubjectName); 
disp(['SubjectName = ', SubjectName])
% Sets the paths and directories
subjpath = pwd; 
load t2stack.mat
load t2mask_final.mat
%load imgraw_ep2d.mat
!cp -r NLI_Outputs/ NLI_OutputsA1
!cp mre_for_inversion.mat mre_for_inversionA1.mat

cd File_Storage/
dir1 = dir('Ax_BRAIN_MRE/*.dcm');
info = dicominfo(fullfile('Ax_BRAIN_MRE', dir1(1).name));
cd ./
load maskx.mat
save maskxA1.mat maskx
% Section 1
dx = (info.ReconstructionDiameter)/double(info.Height);
dy = (info.ReconstructionDiameter)/double(info.Width);
dz = info.SliceThickness; 
freq =double(50); % TS: We need to check where the frequency is in the header

cd('Ax_BRAIN_MRE')
[phsimg] = GE_MRE_RemaskProcessingPart1(dx,dy,dz,freq); 

%% Manual Masking (Human Brain)
close all;
%load('t2mask_bet.mat')
%load('t2stack.mat')
tmp = (t2stack.*mask)./(8000);

% look at the whole brain
figure;im(tmp);caxis([0 1])
while true
    % Display the figure to allow user to view the mask and slices
    figure;
    im(tmp .* mask .* abs(1 - maskx)); 
    caxis([0 1]);

    % Prompt user for the start and end slice
    start_slice = input('Enter the start slice: ');
    end_slice = input('Enter the end slice: ');

    disp(['Start slice: ', num2str(start_slice)]);
    disp(['End slice: ', num2str(end_slice)]);

    last_ss = start_slice;  % Initialize last_ss for error handling

    % Loop through slices in reverse order
    for ss = last_ss:-1:end_slice
        try
            disp(['Processing slice: ', num2str(ss)]);
            
            % Process the mask for the current slice
            maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
            
        catch ME
            % Handle the error, restart the loop from the failed slice
            disp(['Error occurred at slice ', num2str(ss)]);
            disp('Restarting from the last successful slice...');
            
            % Display error message for debugging
            disp(ME.message);
            
            % Set last_ss to the current slice so that the loop restarts from here
            last_ss = ss;
            break;  % Break out of the for loop to restart from the error slice
        end
    end

    % Save the mask after processing all slices
    save('maskx.mat', 'maskx');

    % Display the final mask to check if it's correct
    figure; 
    im(tmp .* mask .* abs(1 - maskx)); 
    caxis([0 1]);

    % Ask the user whether they are satisfied with the mask
    continue_input = input('Please enter Y if you are satisfied with your mask or N if you would like to make changes: ', 's');
    
    % If user is  satisfied, break out of the while loop
    if strcmpi(continue_input, 'Y' )
        disp('User is satisfied with the mask.');
        break;
    % If user wants to make changes, restart the slice selection process
    elseif strcmpi(continue_input, 'N')
        disp('User wants to make changes. Re-enter slices.');
    % Handle invalid inputs
    else
        disp('Invalid input. Please enter Y or N.');
    end
end

%% Section 2
[mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR] = GE_MRE_ProcessingPart2(phsimg,t2stack,mask,dx,dy,dz,freq);

%% Section 3
% clean up data
%LiptonCleanData(mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR)
% Copy To NLI Folder

cd ..
save(sprintf('mre_for_inversion.mat'),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack', 'OSS_SNR')

%addpath(fullfile('/Volumes/McIlvainDrive2/Lipton_Soccer_Study/SUBJECT_DATA', SubjectName))
addpath(fullfile('/Volumes/McIlvainDrive2/Lipton_Lifespan/SUBJECT_DATA', SubjectName))
load('mre_for_inversion.mat')
save(sprintf('%s.mat',SubjectName),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
UIUC_data_convert_mcilvain(SubjectName)
cd(sprintf('%s',SubjectName))
MRE_preprocess_v9_mcilvain('default',SubjectName)
eval(sprintf('!mv %s.mat %s/',SubjectName,SubjectName))
cd ..
%eval(sprintf('!mv %s.mat /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))
%eval(sprintf('!mv %s /Volumes/McIlvainDrive2/Send_to_NLI',SubjectName))

% Ready for Testing (in liu of lines 34 and 35):

%TS 9/2/25 changed to my account because I already had a passkey on the
%computer we need to change it back 
%system(sprintf('scp -r %s gm3128@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName)); 

system(sprintf('scp -r %s ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/',SubjectName)); 
pause(20)
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];
%system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))
system(sprintf('ssh ts3641@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"',insomniapath))

% TS added 11/4/25 for log files and sending pulling back from insomnia23
delete('log/Remasking.txt')

[~, SubjectName] = system('basename "$PWD"');

logFilePath = fullfile('log/Remasked.txt');
fid = fopen(logFilePath, 'a');
if fid ~= -1
    fprintf(fid, 'Remasked NLI for subject: %s\n', SubjectName);
    fclose(fid);
end
logFilePath = fullfile('log/PreNLI_complete.txt');
fid = fopen(logFilePath, 'a');
if fid ~= -1
    fprintf(fid, 'Pre-NLI processing complete for subject: %s\n', SubjectName);
    fclose(fid);
end

%% Section 5
movefile NLI_OutputsA1/ BadMask
movefile mre_for_inversionA1.mat BadMask

movefile BadMask Archive/