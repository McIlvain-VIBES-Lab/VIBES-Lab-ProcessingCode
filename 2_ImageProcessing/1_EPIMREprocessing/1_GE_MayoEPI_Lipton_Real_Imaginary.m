%% Use this code starting with only the Ax_BRAIN_MRE folder
%Section 1

%Comment this back in if you have a T1 Anatomical:
%    cd('T1W_3D_TFE')
%    !/Applications/MRIcron/dcm2nii * -4fn
%    !gunzip -f *.nii.gz
%    !mv co*.nii coT1W_3D_TFE.nii
%    !cp coT1W_3D_TFE.nii ../coT1W_3D_TFE.nii
%    cd ..
%%


%flags for dcm2nii: -f %d %s
cd('Ax_BRAIN_MRE')
!/Applications/MRIcron/dcm2nii  -r y -f %d -n y -o . *
dirlist = dir('*_real.nii');
dirlist2 = dir('*_imaginary.nii');
e1_r = load_nii(sprintf('%s',dirlist(1).name));
e1_i = load_nii(sprintf('%s',dirlist2(1).name));
e1 = double(e1_r.img)+i*double(e1_i.img);

e2_r = load_nii(sprintf('%s',dirlist(2).name));
e2_i = load_nii(sprintf('%s',dirlist2(2).name));
e2 = double(e2_r.img)+i*double(e2_i.img);

e3_r = load_nii(sprintf('%s',dirlist(3).name));
e3_i = load_nii(sprintf('%s',dirlist2(3).name));
e3 = double(e3_r.img)+i*double(e3_i.img);

e4_r = load_nii(sprintf('%s',dirlist(4).name));
e4_i = load_nii(sprintf('%s',dirlist2(4).name));
e4 = double(e4_r.img)+i*double(e4_i.img);

e5_r = load_nii(sprintf('%s',dirlist(5).name));
e5_i = load_nii(sprintf('%s',dirlist2(5).name));
e5 = double(e5_r.img)+i*double(e5_i.img);

e6_r = load_nii(sprintf('%s',dirlist(6).name));
e6_i = load_nii(sprintf('%s',dirlist2(6).name));
e6 = double(e6_r.img)+i*double(e6_i.img);


% 
a = load_nii('18991230_000000s002a1001.nii');
b = load_nii('18991230_000000s002a1001.nii');
a = a.img;

e1_r = a(:,:,:,1:4);
e1_i = a(:,:,:,5:8);
e1 = double(e1_r)+i*double(e1_i);

e2_r = a(:,:,:,9:12);
e2_i = a(:,:,:,13:16);
e2 = double(e2_r)+i*double(e2_i);

e3_r = a(:,:,:,17:20);
e3_i = a(:,:,:,21:24);
e3 = double(e3_r)+i*double(e3_i);

e4_r = a(:,:,:,25:28);
e4_i = a(:,:,:,29:32);
e4 = double(e4_r)+i*double(e4_i);

e5_r = a(:,:,:,33:36);
e5_i = a(:,:,:,37:40);
e5 = double(e5_r)+i*double(e5_i);

e6_r = a(:,:,:,41:44);
e6_i = a(:,:,:,45:48);
e6 = double(e6_r)+i*double(e6_i);

cd ..

%%

dx = b.hdr.dime.pixdim(2);
dy = b.hdr.dime.pixdim(3);
dz = b.hdr.dime.pixdim(4); 

freq = 50;

phsimg(:,:,:,3,:)=angle(e1./e2);
magimg(:,:,:,3,:)=((abs(e1)+abs(e2))/max(col(abs(e1)+abs(e2))))*30000;
phsimg(:,:,:,2,:)=angle(e3./e4);
magimg(:,:,:,2,:)=((abs(e3)+abs(e4))/max(col(abs(e3)+abs(e4))))*30000;
phsimg(:,:,:,1,:)=angle(e5./e6);
magimg(:,:,:,1,:)=((abs(e5)+abs(e6))/max(col(abs(e5)+abs(e6))))*30000;

t2stack = mean(mean(magimg,5),4);
t2nii = make_nii(t2stack,[dy dx dz]);
save_nii(t2nii,'t2stack.nii')

!$FSLDIR/bin/bet2 t2stack.nii t2bet.nii -m -v -f 0.25 -w 1.1
!gunzip -f t2bet.nii_mask.nii.gz
!gunzip -f t2bet.nii.gz
!cp t2bet.nii_mask.nii t2mask.nii
!rm t2bet.nii_mask.nii

tmp = load_nii('t2mask.nii');
seD1 = strel('diamond',0); 
%imerode gets ride of the extra one bad voxel at the edge
mask = imerode(double(permute(flip(flip(tmp.img,1),2),[2 1 3])),seD1);
save t2mask_bet.mat mask

t2stack = flip(flip(permute(t2stack,[2 1 3]),1),2);
save t2stack.mat t2stack


%% Manual Masking (Human Brain)

load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(6000);

maskx = zeros(size(mask));

% look at the whole brain
figure;im(tmp);caxis([0 1])

% look at the mask you have created (the negative mask)
figure;im(tmp.*mask.*abs(1-maskx));caxis([0 1])


for ss = 17:-1:1
    ss
    maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
end
save maskx.mat maskx 

%% OR for Phantoms

maskx = zeros(size(mask));
save maskx.mat maskx 

%% Section 3

load maskx.mat

mask = mask.*abs(1-maskx);
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

save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')


% mkdir(sprintf('%s_Caviness',date))
% cd(sprintf('%s_Caviness',date))
% save(sprintf('%s.mat',mreParams.subj),'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR')
% cd ..

% Clean Up data (Section 4)
mkdir('File_Storage')
% !mv FE_Phs_and_Mag/ File_Storage
% !mv PE_Phs_and_Mag/ File_Storage
% !mv SS_Phs_and_Mag/ File_Storage
!mv Ax_BRAIN_MRE/ File_Storage/
!mv study/ File_Storage
!mv T1W_3D_TFE/ File_Storage/
!mv dcfiles/ File_Storage/

!mv maskx.mat File_Storage/
!mv mre_mag.nii File_Storage/
!mv mre_phs.nii File_Storage/
!mv mre_mask.nii File_Storage/
!mv mre_output.nii File_Storage/

!mv mreimages_unwrap.mat File_Storage/
!mv t2mask_bet.mat File_Storage/