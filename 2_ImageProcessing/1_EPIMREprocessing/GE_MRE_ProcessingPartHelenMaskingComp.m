%function [mreParams,mask,Zmotion,Ymotion,Xmotion,t2stack,OSS_SNR]= GE_MRE_ProcessingPart2


cd File_Storage/Ax_BRAIN_MRE/
%cd Ax_BRAIN_MRE/
dirlist2 = dir('IM*');
nn = 0;
for bb =1:4
for aa = 1:(length(dirlist2)/4) %% number of slices
%    for dd = 1:45
        nn=nn+1;
        images1(:,:,aa,bb) = dicomread(dirlist2(nn).name);
    end
end

tt = 1;
rr = 1;
for jj=1:(length(dirlist2)/48)
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

phsimg(:,:,:,3,:)=angle(e1./e2);
magimg(:,:,:,3,:)=((abs(e1)+abs(e2))/max(col(abs(e1)+abs(e2))))*30000;
phsimg(:,:,:,2,:)=angle(e3./e4);
magimg(:,:,:,2,:)=((abs(e3)+abs(e4))/max(col(abs(e3)+abs(e4))))*30000;
phsimg(:,:,:,1,:)=angle(e5./e6);
magimg(:,:,:,1,:)=((abs(e5)+abs(e6))/max(col(abs(e5)+abs(e6))))*30000;

cd ..

%Section 2 of the MRE Processing Code
% GE_MRE_Processing Part 2
% March 26th 2025
% Grace McIlvain
dx = 1.875;
dy = 1.875;
dz = 2.5;
freq = 50;
% cd ..
load File_Storage/t2mask_bet.mat
load t2stack.mat
%load File_Storage/mreimages_unwrap.mat

% This code makes phase into wave maps for NLI
%load maskx.mat

%%
load('G014_ZS_mask.mat')
mask = mask.*abs(1-maskx);
save t2mask_final.mat mask
% mask = load('t2mask_HL_final.mat');
% mask = mask.mask;
% save t2mask_final.mat mask
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
SubjectName ='Helen5-2023-U7487-0706-VP-01_ZS';


% Save variables to a new .mat for this mask version
save(sprintf('%s.mat', SubjectName), 'mreParams','mask','Zmotion','Ymotion','Xmotion','t2stack','OSS_SNR');

% Convert to NLI-compatible data format
UIUC_data_convert_mcilvain(SubjectName);

% Run preprocessing
cd(SubjectName);
MRE_preprocess_v9_mcilvain('default', SubjectName);
eval(sprintf('!mv %s.mat %s/', SubjectName, SubjectName));
cd ..

% Transfer to NLI (Grace/Teah/Hailin path setup)
system(sprintf('scp -r %s ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName));
pause(20);

% Define path for submissio
insomniapath = ['/insomnia001/depts/mcilvain/users/mcilvain/', SubjectName, '/hex/', SubjectName, '_voxelmesh'];

% % Submit to SLURM via both users (in case)
% system(sprintf('ssh gm3128@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"', insomniapath));
system(sprintf('ssh ts3641@insomnia.rcs.columbia.edu "cd /%s/ && sbatch McIlvain-Submitv9_visc_incomp"', insomniapath));
% 

%end 


%% Check for NLI Outputs
code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/LiptonPipeline';
common_code_path = '/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/2_ImageProcessing/1_EPIMREprocessing/';
addpath('/Volumes/McIlvainDrive2/VIBES-Lab-ProcessingCode/1_StartupCode');
addpath(code_path);
addpath(common_code_path)
startup_matlab_general


SubjectName ='Helen5-2023-U7487-0706-VP-01_ZS';
system(sprintf('scp -r ts3641@insomnia.rcs.columbia.edu:/insomnia001/depts/mcilvain/users/mcilvain/%s .', SubjectName));
cd(SubjectName)
MRE_v9_process_folder

cd hex;
dir2 = dir('*voxelmesh');
cd(dir2.name);
cd inv;
dir2 = dir('*avlast08..0100.Re*');
load(dir2(1).name);
cd ../../../
Rx = double(abs(RealShear-3300)<0.0001);
Rxm = abs(1-Rx);
RealShear = RealShear.*Rxm;
ImagShear = ImagShear.*Rxm;
DR = DR.*Rxm;
save RealShear.mat RealShear
save ImagShear.mat ImagShear
save DR.mat DR

ComplexShear = RealShear + (i*ImagShear);
AbsShear = sqrt(((RealShear.^2)+(ImagShear.^2)));
Mu = 2*(AbsShear.^2)./(RealShear+AbsShear);
save ComplexShear.mat ComplexShear
save AbsShear.mat AbsShear
save Mu.mat Mu
figure;im(Mu(:,:,:)); caxis([0 6000]); colorbar; colormap(gca,stiff_color);
print('-dpng','-r300',sprintf('Mu_%s',dirlist(ii).name(1:end)))