%% Run the starter code

clear all; close all;
delete(gcp('nocreate'));
disp('Code Starting Please Wait')

code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = ['/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/'];
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general

%%
% cd('FE_Phs_and_Mag')
%  !/Applications/MRIcron/dcm2nii * -4fn
%             !gunzip -f *.nii.gz
%             !mv *.nii data.nii
%             cd ..
% cd('PE_Phs_and_Mag')
%  !/Applications/MRIcron/dcm2nii * -4fn
%             !gunzip -f *.nii.gz
%             !mv *.nii data.nii
%             cd ..
% cd('SS_Phs_and_Mag')
%  !/Applications/MRIcron/dcm2nii * -4fn
%             !gunzip -f *.nii.gz
%             !mv *.nii data.nii
%             cd ..
%% Masking

% load('t2mask_bet.mat')
% load('t2stack.mat')
% tmp = (t2stack.*mask)./(10000);
% 
% maskx = zeros(size(mask));
% 
% % look at the whole brain
% figure;im(tmp);caxis([0 1])
% 
% % look at the mask you have created (the negative mask)
% figure;im(tmp.*mask.*abs(1-maskx));caxis([0 1])
% 
% 
% for ss = 20:-1:14
%     ss
%     maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
% end
% save maskx.mat maskx 

%% make inversion files 

clear all
tx = load_nii('FE_Phs_and_Mag/data.nii');
ty = load_nii('PE_Phs_and_Mag/data.nii');
tz = load_nii('SS_Phs_and_Mag/data.nii');

rdx(:,:,:,:,3) = tx.img;
rdx(:,:,:,:,2) = ty.img;
rdx(:,:,:,:,1) = tz.img;

dx = tx.hdr.dime.pixdim(2);
dy = tx.hdr.dime.pixdim(3);
dz = tx.hdr.dime.pixdim(4)/2; % divide by 2 because of the weird way dcm2nii works here

freq = 50;
% currently hardcoded for their 39 slices, not sure how the logic would
% have to be adjusted for other slice numbers
for jj = 1:3
for ii = 1:4 %phase offsets
    
    % phsimg(:,:,1:2:48,jj,ii) = double(rdx(:,:,1:24,ii,jj))./1000;
    % phsimg(:,:,2:2:48,jj,ii) = double(rdx(:,:,1:24,ii+4,jj))./1000;
    % 
    % magimg(:,:,2:2:48,jj,ii) = double(rdx(:,:,25:48,ii,jj));
    % magimg(:,:,1:2:48,jj,ii) = double(rdx(:,:,25:48,ii+4,jj));
% %     
     magimg(:,:,2:2:48,jj,ii) = double(rdx(:,:,1:24,ii,jj));
     magimg(:,:,1:2:48,jj,ii) = double(rdx(:,:,1:24,ii+4,jj));

     phsimg(:,:,2:2:48,jj,ii) = double(rdx(:,:,25:48,ii,jj))./1000;
     phsimg(:,:,1:2:48,jj,ii) = double(rdx(:,:,25:48,ii+4,jj))./1000;
end
end
save magimg magimg

t2stack = mean(mean(magimg,5),4);
t2nii = make_nii(t2stack,[dy dx dz]);
save_nii(t2nii,'t2stack.nii')

!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 1.1
!gunzip -f t2bet.nii_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet.nii_mask.nii t2mask.nii
!rm t2bet.nii_mask.nii

%%

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',0); 
%imerode gets ride of the extra one bad voxel at the edge
mask = imerode(double(permute(flip(flip(tmp.img,1),2),[2 1 3])),seD1);
save t2mask_bet.mat mask

t2stack = flip(flip(permute(t2stack,[2 1 3]),1),2);
save t2stack.mat t2stack

%% 
mask = load('prediction_sample_05_G044.mat');
SubjectName ='prediction_sample_00_G044';
mask = mask.data;
mask = squeeze(mask);

% mask = mask.data;
% squeeze mask
% mask = squeeze(mask);
% mask = mask.*abs(1-maskx); % use this if the mask created is the inverse 
save t2mask_final.mat mask


%load maskx.mat

%mask = mask.*abs(1-maskx);
save t2mask_final.mat mask

image = flip(flip(permute(phsimg,[2 1 3 4 5]),2),1);


 image(:,:,:,2,:) = (-1)*image(:,:,:,2,:);
%image = (-1)*image;

[ny,nx,nz,ni,nj] = size(image);
nk = ni*nj;
image(isnan(image))=0;
imrs = reshape(image,[ny nx nz nk]);

mask_rep = repmat(mask,[1 1 1 nk]);
mask_nii = make_nii(mask_rep);
save_nii(mask_nii,'mre_mask.nii')

phs_nii = make_nii(mask_rep.*imrs);
save_nii(phs_nii,'mre_phs.nii')

mag_rep = repmat(t2stack,[1 1 1 nk]);
mag_nii = make_nii(mag_rep);
save_nii(mag_nii,'mre_mag.nii')

disp('running prelude in fsl... ')
tic
!$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii -v
!gunzip -f mre_output.nii.gz
toc

tmp = load_nii('mre_output.nii');
pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
save mreimages_unwrap.mat pimg
OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);

%disp_img = flip(pimg,5)*1.576; %try this to remove hotness maybe
disp_img = pimg*1.576; % scaling for um, unimportant
fft_img = fft(disp_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

Zmotion = wave_img(:,:,:,1);
Ymotion = wave_img(:,:,:,2);
Xmotion = wave_img(:,:,:,3);

mreParams.subj = 'mre_for_inversion';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = ny;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
