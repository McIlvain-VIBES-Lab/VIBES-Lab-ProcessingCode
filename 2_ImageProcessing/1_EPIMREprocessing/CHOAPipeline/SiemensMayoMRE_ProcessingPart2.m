function [mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR]= SiemensMayoMRE_ProcessingPart2(phsimg,t2stack,mask,dx,dy,dz,freq)

%% Section 2 of the Siemens Mayo MRE Processing Code
% March 26th 2025
% Grace McIlvain

load('imgraw_ep2d.mat') 

mkdir mre_process
cd mre_process

%image(:,:,:,1,:) = squeeze(imgraw(:,:,:,3,:,1)./imgraw(:,:,:,3,:,2)); % Z
%image(:,:,:,2,:) = squeeze(imgraw(:,:,:,2,:,1)./imgraw(:,:,:,2,:,2)); % Y
%image(:,:,:,3,:) = squeeze(imgraw(:,:,:,1,:,1)./imgraw(:,:,:,1,:,2)); % X
image=squeeze(imgraw);
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

end 