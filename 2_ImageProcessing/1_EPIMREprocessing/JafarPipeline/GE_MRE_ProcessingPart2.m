function [mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR]= GE_MRE_ProcessingPart2(phsimg,t2stack,mask,dx,dy,dz,freq)

%% Section 2 of the MRE Processing Code
% GE_MRE_Processing Part 2
% March 26th 2025
% Grace McIlvain

% This code makes phase into wave maps for NLI
load('t2mask_bet.mat')

if exist('maskx.mat')
load maskx.mat
mask = mask.*abs(1-maskx);
end
save t2mask_final.mat mask

image = flip(flip(permute(phsimg,[2 1 3 4 5]),2),1);
image(:,:,:,2,:) = (-1)*image(:,:,:,2,:);

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
!cp dcfiles/OSS_SNR.mat ./

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

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack', 'OSS_SNR')


% mkdir(sprintf('%s_Caviness',date))
% cd(sprintf('%s_Caviness',date))
% save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
cd ..

end 