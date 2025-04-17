function ep2d_v2

delete(gcp('nocreate')); % Clos eparalell pools
% hardwired MRE params
dirs = 3;
posneg =2;
offset = 8;
freq = 25%48.5; % CHANGE THIS!!!!
dirtyp=1; % 1 for curtis (dont change slice orientation or phase encode direction), 2 for matt (try to account for slice orientation/ phase encode)

dirlist_mag = dir('mag/*.IMA');
dirlist_phs = dir('phs/*.IMA');
dcmhead = dicominfo(sprintf('mag/%s',dirlist_mag(1).name));

ny=dcmhead.Rows;
nx=dcmhead.Columns;
nz=length(dirlist_mag)/dirs/posneg/offset; % hardwired for now
dx=dcmhead.PixelSpacing(1);
dy=dcmhead.PixelSpacing(2);
dz=dcmhead.SliceThickness;



totalimg = nz*dirs*posneg*offset;

maglist=dir('mag/*IMA');
phslist=dir('phs/*IMA');

% for ii = 1:totalimg
%     tmp = dir(sprintf('mag/*CURTIS.%04d*.IMA',ii));
%     mag = dicomread(sprintf('mag/%s',tmp.name));
%     tmp = dir(sprintf('phs/*CURTIS.%04d*.IMA',ii));
%     phs = dicomread(sprintf('phs/%s',tmp.name));
%     
%     mag = double(mag);
%     phs = pi*((double(phs)-2048)/2048);
%     cplx_img(:,:,ii) = mag.*exp(1i*phs);
% end

for ii = 1:totalimg
    %tmp = dir(sprintf('mag/*CURTIS.%04d*.IMA',ii));
    mag = dicomread(sprintf('mag/%s',maglist(ii).name));
    %tmp = dir(sprintf('phs/*CURTIS.%04d*.IMA',ii));
    phs = dicomread(sprintf('phs/%s',phslist(ii).name));
    
    mag = double(mag);
    phs = pi*((double(phs)-2048)/2048);
    cplx_img(:,:,ii) = mag.*exp(1i*phs);
end

imgraw = reshape(cplx_img,[ny nx nz posneg dirs offset]);

save imgraw_ep2d.mat imgraw

t2stack = mean(mean(mean(abs(imgraw),6),5),4);
t2nii = make_nii(flipdim(flipdim(permute(t2stack,[2 1 3]),1),2),[dy dx dz]);
save_nii(t2nii,'t2stack.nii')

!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.3 -w 1.3
!gunzip -f t2bet.nii_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet.nii_mask.nii t2mask.nii
!rm t2bet.nii_mask.nii

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',2);
mask = imerode(double(permute(flipdim(flipdim(tmp.img,1),2),[2 1 3])),seD1);
save t2mask_bet.mat mask
save t2stack.mat t2stack

mkdir mre_process
cd mre_process

if(dirtyp==1)
    % Axial scan params from curtis. This goes with a hard-coded Motion2image
    % from UIUC_Data_Convert/MRE_MotionDate_Convert of Motion2Image =
    %      0     1     0
    %     -1     0     0
    %      0     0     1
    image(:,:,:,1,:) = squeeze(imgraw(:,:,:,1,3,:)./imgraw(:,:,:,2,3,:)); % Z
    image(:,:,:,2,:) = (-1)*squeeze(imgraw(:,:,:,1,1,:)./imgraw(:,:,:,2,1,:)); % Y
    image(:,:,:,3,:) = squeeze(imgraw(:,:,:,1,2,:)./imgraw(:,:,:,2,2,:)); % X
elseif(dirtyp==2)
    % Just do nothing and account for in Motion2Image
    image(:,:,:,1,:) = squeeze(imgraw(:,:,:,1,1,:)./imgraw(:,:,:,2,1,:)); % Z
    image(:,:,:,2,:) = squeeze(imgraw(:,:,:,1,2,:)./imgraw(:,:,:,2,2,:)); % Y
    image(:,:,:,3,:) = squeeze(imgraw(:,:,:,1,3,:)./imgraw(:,:,:,2,3,:)); % X
end
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

!$FSLDIR/bin/prelude -a mre_mag.nii -p mre_phs.nii -o mre_output.nii -m mre_mask.nii
!gunzip -f mre_output.nii.gz

tmp = load_nii('mre_output.nii');
pimg = double(reshape(tmp.img,[ny nx nz ni nj]));
save mreimages_unwrap.mat pimg

OSS_SNR = oss_snr_filter(pimg,[dx dy dz],freq,mask);

disp_img = flip(pimg,5)*1.764;
fft_img = fft(disp_img,[],5)/(nj/2);
wave_img = fft_img(:,:,:,:,2);

if(dirtyp==1) 
    Zmotion = wave_img(:,:,:,1);
    Ymotion = wave_img(:,:,:,2);
    Xmotion = wave_img(:,:,:,3);
elseif(dirtyp==2)
    Xmotion = wave_img(:,:,:,1);
    Ymotion = wave_img(:,:,:,2);
    Zmotion = wave_img(:,:,:,3);
end
mreParams.subj = 'mre_for_inversion';
mreParams.FOVx = nx*dx;
mreParams.FOVy = ny*dy;
mreParams.FOVz = nz*dz;
mreParams.nx = nx;
mreParams.ny = ny;
mreParams.nz = nz;
mreParams.freq = freq;
mreParams.oss_snr = OSS_SNR;

if(dirtyp==1)
    save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')
    cd ../
    
    save t2mask.mat mask
    save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack')
elseif(dirtyp==2) % COuldnt get this to reliably work. Sticking with curtis hard-coded version. Wont work if you change the phase encode or slice orientation. But you shouldnt need to. 
    cd ../
    Ur=zeros([size(Xmotion) 3]);
    Ur(:,:,:,1)=real(Xmotion);
    Ur(:,:,:,2)=real(Ymotion);
    Ur(:,:,:,3)=real(Zmotion);
    Ui=zeros([size(Xmotion) 3]);
    Ui(:,:,:,1)=imag(Xmotion);
    Ui(:,:,:,2)=imag(Ymotion);
    Ui(:,:,:,3)=imag(Zmotion);
    [Motion2Image]= Motion2ImageCalc(dcmhead,2);
    AnatomicalMRI=t2stack;
    voxsize_mm=[dx dy dz];
    freqHz=freq;
    RHcoord=true;
    save NLI_MRE_Input Ur Ui AnatomicalMRI RHcoord freqHz voxsize_mm mreParams Motion2Image
    save Mask mask 
end
    
