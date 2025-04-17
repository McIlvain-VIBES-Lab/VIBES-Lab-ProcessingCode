function [images_phs,t2stack] = SiemensMayoMRE_ProcessingPart1(dx,dy,dz,freq)
%% Section 1 of the Siemens Mayo MRE Processing Code
% March 26th 2025
% Grace McIlvain

% This code pulls all dicoms and makes t2stack

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
for dd = 1:24
        for ss = 1:nslice
        sx = slist(ss);
                nn=nn+1;
            images_mag(:,:,sx,dd) = dicomread(dirlist2(nn).name);
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
for dd = 1:24
        for ss = 1:nslice
        sx = slist(ss);
                nn=nn+1;
            images_phs(:,:,sx,dd) = dicomread(dirlist2(nn).name);
        end
    
end
cd ..

images_mag = permute(images_mag,[1 2 3 5 4]);
images_phs = permute(images_phs,[1 2 3 5 4]);

[ny, nx, nz] = size(images_mag(:,:,:,1,1));
% dy = 200/ny;
% dx = 200/nx;
% dz = dy;

% hardwired MRE params (how can we get them to write into the header?) 
dirs = 3; %number of directions (y, x, z) 
posneg = 1; %positive or negative polarity 
offset = 8; %phase offsets 
freq = 60; 

mag2 = double(images_mag);
phs2 = pi*((double(images_phs)-2048)/2048); 
%complex signal is magnitude*exp of the phase 
cplx_img = mag2.*exp(1i*phs2);

imgraw = reshape(cplx_img,[ny nx nz posneg offset dirs]); %offset is the timepoints
imgraw = permute(imgraw,[1 2 3 4 6 5]);
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


end