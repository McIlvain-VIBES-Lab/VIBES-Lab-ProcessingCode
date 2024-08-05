setenv( 'FSLDIR', '/usr/local/fsl' ); % Setting environment variable
setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ' ); % Setting environment variable
fsldir = getenv('FSLDIR'); % gets environment of '/usr/local/fsl', since we have assigned FSLDIR to this value, and assigns it to variable 'fsldir'
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath); % set new path 'fsldirmpath'
setenv( 'ANTSDIR', '/opt/ANTs/bin' ); % Setting environment variable
antsdir = getenv('ANTSDIR');
antsdirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, antsdirmpath); % set new path 'fsldirmpath'
clear fsldir fsldirmpath;

addpath(sprintf('%s/General_Code',pwd))

addpath(sprintf('%s/General_Code/recon_3Dspiral',pwd))
addpath(sprintf('%s/General_Code/recon_3Dspiral/Gridding_for3D',pwd))
addpath(sprintf('%s/General_Code/recon_3Dspiral/nufft',pwd))
addpath(sprintf('%s/General_Code/recon_3Dspiral/irt',pwd))
setup_fessler_IRT(sprintf('%s/General_Code/recon_3Dspiral/irt',pwd))

addpath(genpath(sprintf('%s/General_Code/nifti',pwd)))
addpath(genpath(sprintf('%s/General_Code/oss_snr',pwd)))
cd General_Code
load('mre_colormaps.mat');

cd ..
