%% PROCESS CHOA DATA

% cd ('MRE_3D_AX_ON_AXIS_Mag');
%     !/Applications/MRIcron/dcm2nii * -fn
%     !gunzip -f *.nii.gz
%     !mv *.nii mag.nii
%     %brain extraction tool
%     !$FSLDIR/bin/bet2 mag.nii mag_brain.nii -m -v -f 0.25 -w 1.3
%     !gunzip -f *.nii.gz
%     cd ..
% 
% cd ('MRE_3D_AX_ON_AXIS_P_P');
%     !/Applications/MRIcron/dcm2nii * -fn
%     !gunzip -f *.nii.gz
%     !mv *.nii phs.nii
%     cd .. 

cd('T1_MPRAGE_SAG')
    !/Applications/MRIcron/dcm2nii * -4fn
    !gunzip -f *.nii.gz
    !mv co*.nii coT1_MPRAGE_SAG.nii
    !cp coT1_MPRAGE_SAG.nii ../coT1_MPRAGE_SAG.nii
    cd ..

%%
clear all

cd('MRE_3D_AX_ON_AXIS_Mag')
dirlist2 = dir('IM*');
nimg = length(dirlist2);
nslice = nimg/24; % assumes 24 total images = 8 offsets * 3 directions

% slices are interleaved and not sorted
% assuming siemens convention that starts with 2 for even, and 1 for odd
if ~mod(nslice,2)
    slist = cat(2,2:2:nslice,1:2:nslice);
else
    slist = cat(2,1:2:nslice,2:2:nslice);
end

nn = 0;
for dd = 1:3
    for ss = 1:nslice
        sx = slist(ss);
        for kk = 1:8
            nn=nn+1;
            images_mag(:,:,sx,kk,dd) = dicomread(dirlist2(nn).name);
        end
    end
end

cd ..
cd('MRE_3D_AX_ON_AXIS_P_P')
dirlist2 = dir('IM*');
nimg = length(dirlist2);
nslice = nimg/24; % assumes 24 total images = 8 offsets * 3 directions

% slices are interleaved and not sorted
% assuming siemens convention that starts with 2 for even, and 1 for odd
if ~mod(nslice,2)
    slist = cat(2,2:2:nslice,1:2:nslice);
else
    slist = cat(2,1:2:nslice,2:2:nslice);
end

nn = 0;
for dd = 1:3
    for ss = 1:nslice
        sx = slist(ss);
        for kk = 1:8
            nn=nn+1;
            images_phs(:,:,sx,kk,dd) = dicomread(dirlist2(nn).name);
        end
    end
end
cd ..

images_mag = permute(images_mag,[1 2 3 5 4]);
images_phs = permute(images_phs,[1 2 3 5 4]);

[ny, nx, nz] = size(images_mag(:,:,:,1,1));
dy = 200/ny;
dx = 200/nx;
dz = dy;

% hardwired MRE params (how can we get them to write into the header?) 
dirs = 3; %number of directions (y, x, z) 
posneg = 1; %positive or negative polarity 
offset = 8; %phase offsets 
freq = 60; 

mag2 = double(images_mag);
phs2 = pi*((double(images_phs)-2048)/2048); 
%complex signal is magnitude*exp of the phase 
cplx_img = mag2.*exp(1i*phs2);

imgraw = reshape(cplx_img,[ny nx nz dirs offset posneg]); %offset is the timepoints

save imgraw_ep2d.mat imgraw

t2stack = mean(mean(mean(abs(imgraw),6),5),4);  


t2nii = make_nii(flipdim(flipdim(permute(t2stack,[2 1 3]),1),2),[dy dx dz]);  
save_nii(t2nii,'t2stack.nii')


!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.4
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

load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(2000);

maskx = zeros(size(mask));

% look at the whole brain
figure;im(tmp);caxis([0 1])

for ss = 20:-1:1
    ss
    maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
end
save maskx.mat maskx 


% look at the mask you have created (the negative mask)
figure;im(tmp.*mask.*abs(1-maskx));caxis([0 1])

%% Process Part 2 

mkdir mre_process
cd mre_process

%image(:,:,:,1,:) = squeeze(imgraw(:,:,:,3,:,1)./imgraw(:,:,:,3,:,2)); % Z
%image(:,:,:,2,:) = squeeze(imgraw(:,:,:,2,:,1)./imgraw(:,:,:,2,:,2)); % Y
%image(:,:,:,3,:) = squeeze(imgraw(:,:,:,1,:,1)./imgraw(:,:,:,1,:,2)); % X

image=imgraw;
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

disp_img = (pimg)*1.576;%1.764; TIME REVERAL ON WE TOOK OFF THE FLIP
fft_img = fft(disp_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

Zmotion = wave_img(:,:,:,3);
Ymotion = wave_img(:,:,:,2);
Xmotion = wave_img(:,:,:,1);

mreParams.subj = 'ep2d_slch';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = ny;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;
mreParams.oss_snr = OSS_SNR;

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')

cd ..

save t2mask.mat mask
save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')

disp('finished running ep2d') 

%% Clean Up data
mkdir('File_Storage')
!mv MRE_3D_AX_ON_AXIS_Mag/ File_Storage
!mv MRE_3D_AX_ON_AXIS_P_P/ File_Storage
!mv T1_MPRAGE_SAG/ File_Storage
!mv study/ File_Storage
!mv dcfiles/ File_Storage/

!mv maskx.mat File_Storage/
!mv mre_mag.nii File_Storage/
!mv mre_phs.nii File_Storage/
!mv mre_mask.nii File_Storage/
!mv mre_output.nii File_Storage/

!mv mreimages_unwrap.mat File_Storage/
!mv t2mask_bet.mat File_Storage/




