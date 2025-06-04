%% Dead Brains Wrapper Code
% May 6th 2025
% Some idiot 

clear all; close all;

disp('Code Starting Please Wait')
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
startup_matlab_general

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/JafarPipeline';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing');

addpath(code_path);
subjpath = pwd; 
JafarSetup(subjpath);
addpath('*50*') %JAFAR

dirlist2= dir('*GraceTest') %JAFAR
for jj=1:length(dirlist2)
    display('hi')
    cd(dirlist2(jj).name)
    %dir1 = dir('*.dcm');
    dir1 = dir('dicoms/*.dcm');
    info = dicominfo(fullfile('dicoms', dir1(1).name));
    % Section 1
    dx = (info.ReconstructionDiameter)/double(info.Height);
    dy = (info.ReconstructionDiameter)/double(info.Width);
    dz = info.SliceThickness; 
    freq =200; % JAFAR 
    
    cd('dicoms')
    [phsimg,t2stack] = GE_MRE_ProcessingPart1(dx,dy,dz,freq); 
    
    % Anatomical Scan DCM to NII
    % cd('T1W_3D_TFE')
    % T1_dcm2nii 
    % cd ..
    %% Section 2
    load('t2mask_bet.mat')
    [mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR] = GE_MRE_ProcessingPart2(phsimg,t2stack,mask,dx,dy,dz,freq);
    
    %% Section 3
    cd(dirlist2(jj).name)
    JafarCleanData(mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR)
    cd ..
end