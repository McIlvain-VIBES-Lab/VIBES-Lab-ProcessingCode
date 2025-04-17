% setenv( 'FSLDIR', '/Users/student/fsl' ); % Setting environment variable allen change
setenv( 'FSLDIR', '/Users/vibes-macmini-3/fsl' ); % Setting environment variable
setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ' ); % Setting environment variable
fsldir = getenv('FSLDIR'); % gets environment of '/usr/local/fsl', since we have assigned FSLDIR to this value, and assigns it to variable 'fsldir'
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath); % set new path 'fsldirmpath'
setenv( 'ANTSDIR', '/opt/ANTs/bin' ); % Setting environment variable
antsdir = getenv('ANTSDIR');
antsdirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, antsdirmpath); % set new path 'fsldirmpath'
clear fsldir fsldirmpath;

% TS: added lines 14-24 to run startup code from anywhere
%script_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode';
script_path = mfilename('fullpath');
script_path = script_path(1:end-23);

addpath(sprintf('%s/General_Code',script_path))
clc
addpath(sprintf('%s/General_Code/recon_3Dspiral',script_path))
addpath(sprintf('%s/General_Code/recon_3Dspiral/Gridding_for3D',script_path))
addpath(sprintf('%s/General_Code/recon_3Dspiral/nufft',script_path))
addpath(sprintf('%s/General_Code/recon_3Dspiral/irt',script_path))
setup_fessler_IRT(sprintf('%s/General_Code/recon_3Dspiral/irt',script_path))

addpath(genpath(sprintf('%s/General_Code/nifti',script_path)))
addpath(genpath(sprintf('%s/General_Code/oss_snr',script_path)))

addpath(genpath(sprintf('%s/2_ImageProcessing/4_Inversion/NLI/NLIPrep/MREpreprocess',script_path(1:end-14))))
addpath(genpath(sprintf('%s/2_ImageProcessing/4_Inversion/NLI/NLIPrep/MREpostprocess',script_path(1:end-14))))
addpath(genpath(sprintf('%s/2_ImageProcessing/4_Inversion/NLI/NLIPrep/MREpostprocess/convcode',script_path(1:end-14))))
addpath(genpath(sprintf('%s/2_ImageProcessing/4_Inversion/NLI/NLIPrep/misc',script_path(1:end-14))))

% addpath(sprintf('%s/General_Code',pwd))
% clc
% addpath(sprintf('%s/General_Code/recon_3Dspiral',pwd))
% addpath(sprintf('%s/General_Code/recon_3Dspiral/Gridding_for3D',pwd))
% addpath(sprintf('%s/General_Code/recon_3Dspiral/nufft',pwd))
% addpath(sprintf('%s/General_Code/recon_3Dspiral/irt',pwd))
% setup_fessler_IRT(sprintf('%s/General_Code/recon_3Dspiral/irt',pwd))
% 
% addpath(genpath(sprintf('%s/General_Code/nifti',pwd)))
% addpath(genpath(sprintf('%s/General_Code/oss_snr',pwd)))
cd General_Code
load('mre_colormaps.mat');

cd ..
