function FM = create_field_map(recoInfoFM1)
% function FM = create_field_map(recoInfoFM1)
%
% create field map using Fessler IRT code
% from BPS recon pkg, modified by CLJ
% for use with CLJ "senmap_fm" scans in VB17
%
% NOTE: DELTA_TE is 1 ms --> change here if modified in sequence
%
% needs to have IRT package setup


curdir = pwd;

DELTA_TE = 0.001;

nsl = recoInfoFM1.nsl;
TE1 = recoInfoFM1.TE;
TE2 = TE1 + DELTA_TE;

FM_tp_ind = 1;

% Load in images
if ~exist(sprintf('recon_sen/sen_%05d',FM_tp_ind),'file')
   load(sprintf('recon_sen/sen_%05d',FM_tp_ind))
   Image1 = img;
else
   sprintf('Did not find recon_sen/sen_%05d \n',FM_tp_ind)
   return
end

if ~exist(sprintf('recon_sen/sen_%05d',FM_tp_ind+1),'file')
   load(sprintf('recon_sen/sen_%05d',FM_tp_ind+1))
   Image2 = img;
else
   sprintf('Did not find recon_sen/sen_%05d \n',FM_tp_ind+1)
   return
end


 FM = -mri_field_map_reg(cat(4,Image1, Image2), [TE1 TE2],'l2b',-1,'niter',1000);

cd(curdir)
save FM FM

