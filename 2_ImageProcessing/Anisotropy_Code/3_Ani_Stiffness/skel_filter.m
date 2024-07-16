  function [Mus1,Muf1,Mus2,Muf2,Gp,Gdp] = skel_filter(x,y,z,V1index,V2index,mask,t2stack,dtiall)
tic
[ny, nx, nz, ~] = size(V1index);
ni = 3;
cd 300filt

filt1 = V1index(x,y,z); filt2 = V2index(x,y,z);
tmp1 = load(sprintf('U_filt%d.mat',filt1)); U_filt1 = tmp1.U_filt;
tmp2 = load(sprintf('U_filt%d.mat',filt2)); U_filt2 = tmp2.U_filt;
cd ..

% x1 = U_filt1(:,:,:,1); y1 = U_filt1(:,:,:,2); z1 = U_filt1(:,:,:,3);
% tmpx1 = U_filt1(:,:,:,1)./(x1.^2+y1.^2+z1.^2); % x-component 
% tmpy1 = U_filt1(:,:,:,2)./(x1.^2+y1.^2+z1.^2); % y-component 
% tmpz1 = U_filt1(:,:,:,3)./(x1.^2+y1.^2+z1.^2); % z-component x = 
% 
% x2 = U_filt2(:,:,:,1); y2 = U_filt2(:,:,:,2); z2 = U_filt2(:,:,:,3);
% tmpx2 = U_filt2(:,:,:,1)./(x2.^2+y2.^2+z2.^2);
% tmpy2 = U_filt2(:,:,:,2)./(x2.^2+y2.^2+z2.^2);
% tmpz2 = U_filt2(:,:,:,3)./(x2.^2+y2.^2+z2.^2);
% 
% Ufilt1 = cat(4,tmpx1,tmpy1,tmpz1);
% Ufilt2 = cat(4,tmpx2,tmpy2,tmpz2);

load bucky300.mat
dti = flip(flip(permute(dtiall,[2 1 3 4]),1),2);
a = squeeze(dti(x,y,z,:));
n1 = allnvec(filt1,:)';
n2 = allnvec(filt2,:)';

m_s1 = cross(n1,a)/norm(cross(n1,a));
m_f1 = cross(n1,m_s1)/norm(cross(n1,m_s1));
m_s2 = cross(n2,a)/norm(cross(n2,a));
m_f2 = cross(n2,m_s2)/norm(cross(n2,m_s2));

m_s1x = repmat(permute(m_s1,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
m_s2x = repmat(permute(m_s2,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
m_f1x = repmat(permute(m_f1,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);
m_f2x = repmat(permute(m_f2,[2 3 4 1]),[size(U_filt1,1) size(U_filt1,2) size(U_filt1,3) 1]);

Us1 = repmat(dot(U_filt1,m_s1x,4),[1 1 1 3]).*m_s1x; Uf1 = repmat(dot(U_filt1,m_f1x,4),[1 1 1 3]).*m_f1x; 
Us2 = repmat(dot(U_filt2,m_s2x,4),[1 1 1 3]).*m_s2x; Uf2 = repmat(dot(U_filt2,m_f2x,4),[1 1 1 3]).*m_f2x; 

% figure;imagesc(abs(squeeze(ffp(Us1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Us2(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf2(:,:,24,:))))); axis square

%% LDI PREPPING
tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Us1,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);

userow = 1:size(PCdata,1);
usecol = 1:size(PCdata,2);
subjID = 'subj2_Us1';
nanmask = mask;
nanmask(nanmask == 0) = NaN;
nanmask = nanmask;
voxelsize = 3;
MAG = repmat(t2stack,[1 1 1 4 3]);
PC2micron = 2.6250;
[Gps1,Gdps1] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Us1_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mus1 = abs(Gps1(x,y,z)+1i*Gdps1(x,y,z));
Gps1x = Gps1(x,y,z);
Gdps1x = Gdps1(x,y,z);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Uf1,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf1,Gdpf1] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf1_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Muf1 = abs(Gpf1(x,y,z)+1i*Gdpf1(x,y,z));
Gpf1x = Gpf1(x,y,z);
Gdpf1x = Gdpf1(x,y,z);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Us2,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gps2,Gdps2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Us2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mus2 = abs(Gps2(x,y,z)+1i*Gdps2(x,y,z));
Gps2x = Gps2(x,y,z);
Gdps2x = Gdps2(x,y,z);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(Uf2,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf2,Gdpf2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Muf2 = abs(Gpf2(x,y,z)+1i*Gdpf2(x,y,z));
Gpf2x = Gpf2(x,y,z);
Gdpf2x = Gdpf2(x,y,z);

Gp = [Gps1x Gpf1x Gps2x Gpf2x];
Gdp = [Gdps1x Gdpf1x Gdps2x Gdpf2x];
toc
