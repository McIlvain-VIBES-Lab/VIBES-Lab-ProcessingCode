  function [Mus,Muf] = MuSF(mask,t2stack,dtiall)
tic
ni = 3;
[ny, nx, nz, ~] = size(mask);
load BMR.mat
U3 = Brfo_all + 1i.*Bifo_all;
Ui = zeros(size(U3));
xy = length(U3(:,1,1,1));
z = length(U3(1,1,:,1));
tic
for dir = 1:3
    U3(isnan(U3)) = 0;
    % G3 = spatial FFT of complex coefficient field, needs to be fftshifted to center it
    G3 = fftshift(fftn(U3(:,:,:,dir).*mask));    
    % set up kspace
    [ay,ax,az]=size(G3);    
    % radial and angular coordinates in kspace
    [KX,KY,KZ]=meshgrid(1:ax,1:ay,1:az);
    % center of kspace
    cenx=floor(ax/2)+1;
    ceny=floor(ay/2)+1;
    cenz=floor(az/2)+1;
    KX=KX-cenx;
    KY=KY-ceny;
    KZ=KZ-cenz;        
    % magnitude of wavenumber if any radial filtering
    KR = sqrt(KX.^2+KY.^2+KZ.^2);        
    % % filter in multiple directions spanning spherical angles
    % dth = pi/Ntheta;
    % dph = pi/Ntheta;
    % th00 = [0:dth:(2*pi-dth)];
    % phi00 = [-pi/2:dph:pi/2];    
    % th00 = [0:1:0];
    % phi00 = [0:1:0];
    % load buckyinfo.mat allnvec;
    for nn=1:92;        
        th0 = atan2(allnvec(nn,2),allnvec(nn,1));
        phi0 = atan2(allnvec(nn,3),allnvec(nn,1)/cos(th0));        
        % normal vector
        no0=[cos(th0)*cos(phi0);sin(th0)*cos(phi0);sin(phi0)];        
        % directional cosine
        cn = (no0(1)*KX + no0(2)*KY + no0(3)*KZ)./KR;
        cn(find(isnan(cn)))= zeros(size(find(isnan(cn))));       
        % Define angular filter coefficients
        % cos(4*theta).^2
        % cn4 = (64*cn.^8 - 128*cn.^6 +80*cn.^4 -16*cn.^2+1);
        %cn4 = cn4.*(cn>cos(pi/8));
        % cos(2*theta).^2
        cn2 = (4*cn.^4-4*cn.^2+1);
        cn2 = cn2.*(cn>cos(pi/4));        
        % Develop radial filter here if necessary
        % omit       
        % Apply the filter here
        % filter
        Gf = G3.*cn2;
        %   Uf(:,:,:,nn)=Gf;
        % Uf(:,:,:,nn) = (ifftn(ifftshift(Gf)));
        Ui(:,:,:,dir,nn) = (ifftn(ifftshift(Gf)));
        %%%%%%%%%%%%%%%%%%%%%%%%%  
        % vector fields from weighted amplitudes multiplied by direction
        ux = squeeze(abs(Ui(:,:,:,dir,nn)))*no0(1);
        uy = squeeze(abs(Ui(:,:,:,dir,nn)))*no0(2);
        uz = squeeze(abs(Ui(:,:,:,dir,nn)))*no0(3);
        UXi(:,:,:,dir,nn) = ux.*mask;       % x component
        UYi(:,:,:,dir,nn) = uy.*mask;        % y component
        UZi(:,:,:,dir,nn) = uz.*mask;        % z component
    end  
end  

% weighted average direction of propagation
UUX=-sum(UXi,5); UUY=-sum(UYi,5); UUZ=-sum(UZi,5);
UX = sum(UUX,4); UY = sum(UUY,4); UZ = sum(UUZ,4);
UX_norm = UX./sqrt(UX.^2+UY.^2+UZ.^2); UY_norm = UY./sqrt(UX.^2+UY.^2+UZ.^2); UZ_norm = UZ./sqrt(UX.^2+UY.^2+UZ.^2);

N = cat(4,UX,UY,UZ);
% N = flip(flip(permute(cat(4,UX,UY,UZ),[2 1 3 4]),2),1);

N_norm = cat(4,UX_norm,UY_norm,UZ_norm);
N_norm(N_norm==0) = NaN;

M_s = cross(N_norm,ffp(dtiall),4);
M_sx = M_s(:,:,:,1); M_sy = M_s(:,:,:,2); M_sz = M_s(:,:,:,3);
M_sunit(:,:,:,1) = M_sx./sqrt(M_sx.^2+M_sy.^2+M_sz.^2);
M_sunit(:,:,:,2) = M_sy./sqrt(M_sx.^2+M_sy.^2+M_sz.^2);
M_sunit(:,:,:,3) = M_sz./sqrt(M_sx.^2+M_sy.^2+M_sz.^2);
M_f = cross(N_norm,M_s,4);
M_fx = M_f(:,:,:,1); M_fy = M_f(:,:,:,2); M_fz = M_f(:,:,:,3);
M_funit(:,:,:,1) = M_fx./sqrt(M_fx.^2+M_fy.^2+M_fz.^2);
M_funit(:,:,:,2) = M_fy./sqrt(M_fx.^2+M_fy.^2+M_fz.^2);
M_funit(:,:,:,3) = M_fz./sqrt(M_fx.^2+M_fy.^2+M_fz.^2);

U_sdot = dot(N,M_sunit,4); U_fdot = dot(N,M_funit,4);

for dir = 1:3
    U_s(:,:,:,dir) = M_sunit(:,:,:,dir).*U_sdot;
    U_f(:,:,:,dir) = M_funit(:,:,:,dir).*U_fdot;
end

% figure;imagesc(abs(squeeze(ffp(Us1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf1(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Us2(:,:,24,:))))); axis square
% figure;imagesc(abs(squeeze(ffp(Uf2(:,:,24,:))))); axis square

%% LDI PREPPING
tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(U_s,[1 1 1 1 4]).*exp(1i*tt));
PCdata = permute(ttmp,[1 2 3 5 4]);

userow = 10:size(PCdata,1)-10;
usecol = 10:size(PCdata,2)-10;
subjID = 'subj2_Us1';
nanmask = mask;
nanmask(nanmask == 0) = NaN;
nanmask = nanmask;
voxelsize = 3;
MAG = repmat(t2stack,[1 1 1 4 3]);
PC2micron = 2.6250;

[Gps2,Gdps2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Us2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mus = abs(Gps2+1i*Gdps2);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(U_f,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf2,Gdpf2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Muf = abs(Gpf2+1i*Gdpf2);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(N,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf2,Gdpf2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mu_base = abs(Gpf2+1i*Gdpf2);

tspc = (2*pi/4);
t= (0:tspc:((2*pi)-tspc))';
tt = repmat(reshape(t,[1 1 1 1 4]),[ny nx nz ni 1]);
ttmp = real(repmat(U3,[1 1 1 1 4]).*exp(1i*tt));

PCdata = permute(ttmp,[1 2 3 5 4]);
[Gpf2,Gdpf2] = run_MRE_LDI_script_fnc(PCdata,userow,usecol,mask,nanmask,voxelsize,MAG,PC2micron,subjID);
%save subj2_Uf2_forLDI.mat PCdata userow usecol subjID nanmask voxelsize MAG PC2micron mask
Mu_U3 = abs(Gpf2+1i*Gdpf2);

toc
