% Helen's script for the MRE_Vol_xx data 
% Used to process Mre_Vol_8

cd('mag')
 !/Applications/MRIcron/dcm2nii * -4fn
            !gunzip -f *.nii.gz
            !mv *.nii mag.nii
            cd ..
cd('phs')
 !/Applications/MRIcron/dcm2nii * -4fn
            !gunzip -f *.nii.gz
            !mv *.nii phs.nii
            cd ..


 % cd('T1W_3D_TFE')
 % !/Applications/MRIcron/dcm2nii * -4fn
 %            !gunzip -f *.nii.gz
 %            !mv co*.nii coT1W_3D_TFE.nii
 %            !cp coT1W_3D_TFE.nii ../coT1W_3D_TFE.nii

       %     cd ..
 %% Processing Part 1
%clear all

    addpath(sprintf('%s/mag',pwd))
    addpath(sprintf('%s/phs',pwd))


mag_nii = load_untouch_nii('mag/mag.nii');
phs_nii = load_untouch_nii('phs/phs.nii');
magimg = mag_nii.img;
phsimg = phs_nii.img;

mag = permute((flip(magimg,2)),[2 1 3 4]);
phs = permute((flip(phsimg,2)),[2 1 3 4]);
% figure;im(mag(:,:,:,1));

dirlist_mag = dir('mag/*.nii*');
dirlist_phs = dir('phs/*.nii*');

% nx = length(mag(:,1,1,1));  %number of points
% ny = length(mag(1,:,1,1));
% nz = length(mag(1,1,:,1)); % Can also input manually. nz = 60 for 2mm, 48 for 3mm
[ny,nx,nz,~] = size(mag);
dx = mag_nii.hdr.dime.pixdim(2); %spacing in x (size of each voxel) 
dy = mag_nii.hdr.dime.pixdim(3);
dz = mag_nii.hdr.dime.pixdim(4);

% hardwired MRE params (how can we get them to write into the header?) 
dirs = 3; %number of directions (y, x, z) 
posneg = 2; %positive or negative polarity 
offset = 4; %phase offsets 
freq = 60;  

mag = double(mag);
phs = pi*((double(phs)-2048)/2048); 
%complex signal is magnitude*exp of the phase 
cplx_img = mag.*exp(1i*phs);

imgraw = reshape(cplx_img,[ny nx nz posneg dirs offset]); %offset is the timepoints

save imgraw_ep2d.mat imgraw

t2stack = mean(mean(mean(abs(imgraw),6),5),4);  


t2nii = make_nii(flipdim(flipdim(permute(t2stack,[2 1 3]),1),2),[dy dx dz]);  
save_nii(t2nii,'t2stack.nii')


!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 1.3
!gunzip -f t2bet.nii_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet.nii_mask.nii t2mask.nii
!rm t2bet.nii_mask.nii

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',2); 
mask = imerode(double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3 4])),seD1);
save t2mask_bet.mat mask
save t2stack.mat t2stack

%% Manual Masking
% 
% load('t2mask_bet.mat')
% load('t2stack.mat')
% tmp = (t2stack.*mask)./(1000);
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
% for ss = 14:-1:1
%     ss
%     maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
% end
% save maskx.mat maskx 

%% Process Part 2 
clear image
% load maskx.mat
% Helen's mod:
all_masks = {'S071_AH_mask.mat', 'S071_GM_mask.mat', 'S071_ZS_mask.mat', 'S071_TS_mask.mat'};
t2bet = load('t2mask_bet.mat');

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
        maskx = loaded.maskx;
    
    else
        error('The file %s does not contain a variable named "maskx"', mask_file);
    end

    fprintf('\n=== Processing with manual mask: %s ===\n', mask_label);
   
    mask = t2bet.*abs(1-maskx);
    filename = sprintf('t2mask_%s_final.mat', mask_label);
    save(filename, 'mask'); % it will be referred as "mask" for the following code
    
    
    mkdir mre_process
    cd mre_process
    
    image(:,:,:,:,1) = squeeze(imgraw(:,:,:,1,3,:)./imgraw(:,:,:,2,3,:)); % Z
    image(:,:,:,:,2) = (-1)*squeeze(imgraw(:,:,:,1,1,:)./imgraw(:,:,:,2,1,:)); % Y
    image(:,:,:,:,3) = squeeze(imgraw(:,:,:,1,2,:)./imgraw(:,:,:,2,2,:)); % X
    image = permute(image,[1 2 3 5 4]);
    
    save mreimages.mat image
    
    [ny,nx,nz,ni,nj] = size(image);
    nk = ni*nj;
    image(isnan(image))=0;
    imrs = reshape(image,[ny nx nz nk]);
    
    save t2mask.mat mask
    
    mask_rep = repmat(mask,[1 1 1 nk]);
    mask_nii = make_nii(mask_rep);
    save_nii(mask_nii,'mre_mask.nii')
    
    phs_nii = make_nii(mask_rep.*angle(imrs));
    save_nii(phs_nii,'mre_phs.nii')
    
    mag_rep = repmat(t2stack,[1 1 1 nk]);
    mag_nii = make_nii(mag_rep);
    save_nii(mag_nii,'mre_mag.nii')
    
    
    tic
    !$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii -v
    !gunzip -f mre_output.nii.gz
    toc
    
    tmp = load_nii('mre_output.nii');
    pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
    save mreimages_unwrap.mat pimg
    
    
    OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);
    
    disp_img = flip(pimg,5)*1.764;
    fft_img = fft(disp_img,[],5)/(nj/2);
    wave_img = fft_img(:,:,:,:,2);
    
    Zmotion = wave_img(:,:,:,1);
    Ymotion = wave_img(:,:,:,2);
    Xmotion = wave_img(:,:,:,3);
    
    mreParams.subj = 'ep2d_slch';
    mreParams.FOVx = nx*dx;
    mreParams.FOVy = ny*dy;
    mreParams.FOVz = nz*dz;
    mreParams.nx = nx;
    mreParams.ny = ny;
    mreParams.nz = nz;
    mreParams.freq = freq;
    mreParams.oss_snr = OSS_SNR;
    
    save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
    
    cd ..
    
    save t2mask.mat mask
    save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
    
    disp('finished running ep2d') 


    % Section 3A: Send to NLI (mask version)
    HelenLiptonPrepareForNLI(mreParams, mask, Zmotion, Ymotion, Xmotion, t2stack, OSS_SNR, mask_label);
end


