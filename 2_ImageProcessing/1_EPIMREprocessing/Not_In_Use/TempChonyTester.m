%% Grace's Temp Tester CHONY data sorter

clear all
cd('Ax_BRAIN_MRE')
dirlist2 = dir('IM*');
nn = 0;
for bb =1:4
for aa = 1:540
%    for dd = 1:45
        nn=nn+1;
        images1(:,:,aa,bb) = dicomread(dirlist2(nn).name);
    end
end

tt = 1;
rr = 1;
for jj=1:45
for ii=1:4

e1_r(:,:,rr,ii) = images1(:,:,tt,ii);
e1_i(:,:,rr,ii) = images1(:,:,tt+1,ii);

e2_r(:,:,rr,ii) = images1(:,:,tt+2,ii);
e2_i(:,:,rr,ii) = images1(:,:,tt+3,ii);

e3_r(:,:,rr,ii) = images1(:,:,tt+4,ii);
e3_i(:,:,rr,ii) = images1(:,:,tt+5,ii);

e4_r(:,:,rr,ii) = images1(:,:,tt+6,ii);
e4_i(:,:,rr,ii) = images1(:,:,tt+7,ii);

e5_r(:,:,rr,ii) = images1(:,:,tt+8,ii);
e5_i(:,:,rr,ii) = images1(:,:,tt+9,ii);

e6_r(:,:,rr,ii) = images1(:,:,tt+10,ii);
e6_i(:,:,rr,ii) = images1(:,:,tt+11,ii);

end
tt = tt+12;
rr =rr+1;
end


e1_t = double(e1_r)+i*double(e1_i);
e2_t = double(e2_r)+i*double(e2_i);
e3_t = double(e3_r)+i*double(e3_i);
e4_t = double(e4_r)+i*double(e4_i);
e5_t = double(e5_r)+i*double(e5_i);
e6_t = double(e6_r)+i*double(e6_i);

e1 = flip(flip(permute(e1_t,[2 1 3 4]),2),1);
e2 = flip(flip(permute(e2_t,[2 1 3 4]),2),1);
e3 = flip(flip(permute(e3_t,[2 1 3 4]),2),1);
e4 = flip(flip(permute(e4_t,[2 1 3 4]),2),1);
e5 = flip(flip(permute(e5_t,[2 1 3 4]),2),1);
e6 = flip(flip(permute(e6_t,[2 1 3 4]),2),1);

cd ..

%%

dx = 1.8750;
dy = 1.8750;
dz = 2.500; 

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

%% OR for Phantoms

maskx = zeros(size(mask));
save maskx.mat maskx 
%% OR Manual Masking (Human Brain)

load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(6000);

maskx = zeros(size(mask));

% look at the whole brain
figure;im(tmp);caxis([0 1])

% look at the mask you have created (the negative mask)



for ss = 45:-1:1
    ss
    maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
end
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

%% Clean Up data (Section 4)
mkdir('File_Storage')
!mv FE_Phs_and_Mag/ File_Storage
!mv PE_Phs_and_Mag/ File_Storage
!mv SS_Phs_and_Mag/ File_Storage
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