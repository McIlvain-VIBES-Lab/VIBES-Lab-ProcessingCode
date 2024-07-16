  function [Mus,Muf] = skel_filter2(x,y,z,Index,mask,t2stack,dtiall)
tic
[ny, nx, nz, ~] = size(Index);
ni = 3;
cd 300filt

filt = Index(x,y,z);
tmp = load(sprintf('U_filt%d.mat',filt));
U_filt = tmp.Ufilt;

x = U_filt(:,:,:,1); y = U_filt(:,:,:,2); U_filt(:,:,:,3);
tmpx = U_filt(:,:,:,1)./(x.^2+y.^2+z.^2);
tmpy = U_filt(:,:,:,2)./(x.^2+y.^2+z.^2);
tmpz = U_filt(:,:,:,3)./(x.^2+y.^2+z.^2);

Ufilt = cat(4,tmpx,tmpy,tmpz);
m_s = cross(Ufilt,dtiall,4); m_f = cross(Ufilt,dtiall,4);
Us = repmat(dot(U_filt,m_s,4),[1 1 1 3]).*m_s; Uf = repmat(dot(U_filt,m_f,4),[1 1 1 3]).*m_f; 

% figure;imagesc(abs(squeeze(ffp(Us1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Us2(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf2(:,:,24,:))))); axis square

%% LDI PREPPING
userow = 10:size(PCdata,1)-10;
usecol = 10:size(PCdata,2)-10;
subjID = 'subj2_Us1';
nanmask = mask;
nanmask(nanmask == 0) = NaN;
nanmask = nanmask;
voxelsize = 3;
MAG = repmat(t2stack,[1 1 1 4 3]);
PC2micron = 2.6250;

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Us,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gps2,Gdps2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Us2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mus = abs(Gps2(x,y,z)+1i*Gdps2(x,y,z));

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Uf,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf2,Gdpf2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Muf = abs(Gpf2(x,y,z)+1i*Gdpf2(x,y,z));
toc
