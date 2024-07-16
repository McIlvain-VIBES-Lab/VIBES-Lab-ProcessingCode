function [Ui,U_s_V1,U_f_V1,U_s_V2,U_f_V2]=slow_fast_twopeaks_exp1(Brfo_all,Bifo_all,mask,dtiall,V1index,V2index,A1x,A2x)
% show amplitude weighted direction of wave propagation by directional
% filtering
%  % uses 92 directions based on bucky ball vertices and mid faces
% written by R. Okamoto, Washington University 2014
% inputs:
% U3 - (NxMxP) field of complex coefficients (Re,Im) of fundamental component
% mask - mask array (NxM)
% outputs
% UUX, UUY, UUZ: amplitude-weighted vectors of propagation direction
%% set Parameters and Create k-space Directional filter
% Brfo_all = Output of Bulk Motion Removed Code from WashU (Real Component)
% Bifo_all = Output of Bulk Motion Removed code from WashU (Imaginary
% Component)
A = dtiall;
U3 = Brfo_all + 1i.*Bifo_all;
%U3 = xyz4;
Uo = zeros(size(U3));
xy = length(U3(:,1,1,1));
z = length(U3(1,1,:,1));

U3(isnan(U3)) = 0;
% set up kspace
% [ay,ax,az]=size(U3(:,:,:,1));  
np=size(U3(:,:,:,1)); 
dk=size(U3(:,:,:,1)).^-1;
kmin = dk.*(np/2)*(-1);
kmax = dk.*((np/2)-1);
% radial and angular coordinates in kspace
% [KX,KY,KZ]=meshgrid(1:ax,1:ay,1:az);
% % center of kspace
% cenx=floor(ax/2)+1;
% ceny=floor(ay/2)+1;
% cenz=floor(az/2)+1;
% KX=KX-cenx; KY=KY-ceny; KZ=KZ-cenz;     

[KX,KY,KZ]=meshgrid(kmin(2):dk(2):kmax(2),kmin(1):dk(1):kmax(1),kmin(2):dk(3):kmax(3));
% magnitude of wavenumber if any radial filtering
KR = sqrt(KX.^2+KY.^2+KZ.^2);     

klim = min(abs(kmin)); % extent in k-space for the filter to extend (set to min because smaller dimension has less 
filter = 1./(1+(KR/klim).^(10)); % low pass filtering up to klim
filter(((xy/2)+1),((xy/2)+1),((z/2)+1)) = 0;

tic
G3 = zeros(xy,xy,z,3);
U3filt = zeros(size(G3));
for dir = 1:3
    G3(:,:,:,dir) = fftshift(fftn(U3(:,:,:,dir).*mask));
    G3(:,:,:,dir) = squeeze(G3(:,:,:,dir)).*filter;
    U3filt(:,:,:,dir) = (ifftn(ifftshift(squeeze(G3(:,:,:,dir)))));
end
toc
tic
disp('G3 found')
   

load bucky300.mat allnvec;
%allnvec = X;
Theta0 = zeros(size(allnvec(:,1)));
Phi0 = zeros(size(allnvec(:,1)));
disp('setup complete')
cn2_nn = zeros(xy,xy,z,size(allnvec,1));
tic
for nn=1:length(allnvec);
    if rem(nn,100) == 0
        nn
        toc
    end
	th0 = atan2(allnvec(nn,2),allnvec(nn,1));
	Theta0(nn) = th0;
	phi0 = atan2(allnvec(nn,3),allnvec(nn,1)/cos(th0)); 
	Phi0(nn) = phi0;
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
    cn2_nn(:,:,:,nn) = cn2;
end
tic
for dir = 1:3
    Ui = zeros(xy,xy,z);
    for nn = 1:length(allnvec)
        if rem(nn,100) == 0
            nn
            toc
        end
        Gf = G3(:,:,:,dir).*cn2_nn(:,:,:,nn);
        Ui(:,:,:,nn) = (ifftn(ifftshift(Gf)));
    end
    if dir == 1
        Ux = Ui;
        display('Ux found')
    elseif dir == 2;
        Uy = Ui;
        display('Uy found')
    elseif dir == 3;
        Uz = Ui;
        display('Uz found')
    end
    toc
end
display('U assigned')
toc

U3new = permute(cat(5,Ux,Uy,Uz),[1 2 3 5 4]);
Ui = U3new;

%% Cross products to create m_s & m_f
% N = Uf; %defines the base directional vectors
sprintf('U found')
tic
load bucky300.mat
allnvec_mat = permute(repmat(allnvec,[1 1 xy xy z]),[3 4 5 2 1]);
%load 'allnvec_matrix.mat' allnvec_matrix 
%load in matrix for crossing with A, the fiber direction from DTI
allnvec_mat(abs(allnvec_mat)<=.005) = 0; %matrix isn't zeroed out for some reason so just sets numbers to 0 that should be 0

m_si = zeros(xy,xy,z,3,length(allnvec(:,1)));
m_fi = zeros(size(m_si));

m_sx = zeros(size(m_si(:,:,:,1))); m_sy = zeros(size(m_si(:,:,:,1))); m_sz = zeros(size(m_si(:,:,:,1)));
m_si_unit = zeros(size(m_si));

m_fx = zeros(size(m_si(:,:,:,1))); m_fy = zeros(size(m_si(:,:,:,1))); m_fz = zeros(size(m_si(:,:,:,1)));
m_fi_unit = zeros(size(m_si));

%Individual cross products for each of the 92 directions for m_s and m_f
% m_s = NxA
% m_f = m_sxA
tic
 for nn = 1:length(allnvec(:,1))
     nn
     for ii = 1:xy
         for jj = 1:xy
             for kk = 1:z
                 m_si(ii,jj,kk,:,nn) = cross(squeeze(allnvec_mat(ii,jj,kk,:,nn)),squeeze(A(ii,jj,kk,:)),1); 
                 m_sx(ii,jj,kk,nn) = squeeze(m_si(ii,jj,kk,1,nn)); m_sy(ii,jj,kk,nn) = squeeze(m_si(ii,jj,kk,2,nn)); m_sz(ii,jj,kk,nn) = squeeze(m_si(ii,jj,kk,3,nn)); 
                 m_si_unit(ii,jj,kk,1,nn) = m_sx(ii,jj,kk,nn)./(sqrt(m_sx(ii,jj,kk,nn).^2+m_sy(ii,jj,kk,nn).^2+m_sz(ii,jj,kk,nn).^2));
                 m_si_unit(ii,jj,kk,2,nn) = m_sy(ii,jj,kk,nn)./(sqrt(m_sx(ii,jj,kk,nn).^2+m_sy(ii,jj,kk,nn).^2+m_sz(ii,jj,kk,nn).^2));
                 m_si_unit(ii,jj,kk,3,nn) = m_sz(ii,jj,kk,nn)./(sqrt(m_sx(ii,jj,kk,nn).^2+m_sy(ii,jj,kk,nn).^2+m_sz(ii,jj,kk,nn).^2)); 
                 
                 m_fi(ii,jj,kk,:,nn) = cross(squeeze(allnvec_mat(ii,jj,kk,:,nn)),squeeze(m_si_unit(ii,jj,kk,:,nn)),1);
                 m_fx(ii,jj,kk,nn) = squeeze(m_fi(ii,jj,kk,1,nn)); m_fy(ii,jj,kk,nn) = squeeze(m_fi(ii,jj,kk,2,nn)); m_fz(ii,jj,kk,nn) = squeeze(m_fi(ii,jj,kk,3,nn));
                 m_fi_unit(ii,jj,kk,1,nn) = m_fx(ii,jj,kk,nn)./(sqrt(m_fx(ii,jj,kk,nn).^2+m_fy(ii,jj,kk,nn).^2+m_fz(ii,jj,kk,nn).^2));
                 m_fi_unit(ii,jj,kk,2,nn) = m_fy(ii,jj,kk,nn)./(sqrt(m_fx(ii,jj,kk,nn).^2+m_fy(ii,jj,kk,nn).^2+m_fz(ii,jj,kk,nn).^2));
                 m_fi_unit(ii,jj,kk,3,nn) = m_fz(ii,jj,kk,nn)./(sqrt(m_fx(ii,jj,kk,nn).^2+m_fy(ii,jj,kk,nn).^2+m_fz(ii,jj,kk,nn).^2));
             end
         end
     end
 end
toc
sprintf('Waves found')
tic

%% take dot products and multiply by m_s and m_f to create U vectors
U_s_V1 = zeros(xy,xy,z);
U_s_V2 = zeros(xy,xy,z);
U_f_V1 = zeros(xy,xy,z);
U_f_V2 = zeros(xy,xy,z);
U_si = zeros(size(m_si_unit));
U_fi = zeros(size(U_si)); 
U_s_doti = zeros(size(squeeze(m_si_unit(:,:,:,1,:))));
U_f_doti = zeros(size(U_s_doti));
for nn = 1:length(allnvec(:,1))
    %individual dot products for U_s and U_f
    U_s_doti(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_si_unit(:,:,:,:,nn)),4)); 
    U_f_doti(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_fi_unit(:,:,:,:,nn)),4));
%     %taking the individual dot products and then using it as a scaling
%     %measure for the individual directionality vector
%     for dir = 1:3
%         U_si(:,:,:,dir,nn) = squeeze(m_si_unit(:,:,:,dir,nn)).*U_s_doti(:,:,:,nn);
%         U_fi(:,:,:,dir,nn) = squeeze(m_fi_unit(:,:,:,dir,nn)).*U_f_doti(:,:,:,nn);
%     end
end
toc
sprintf('individual slow and fast wave fields found')

% what is left: pulling out the correct U_s_doti for V1 and V2
for ix=1:xy
    for iy = 1:xy
        for iz = 1:z
            U_s_V1(iy,ix,iz) = squeeze(U_s_doti(ix,iy,iz,V1index(iy,ix,iz)));
            U_f_V1(iy,ix,iz) = squeeze(U_f_doti(ix,iy,iz,V1index(iy,ix,iz)));
            U_s_V2(iy,ix,iz) = squeeze(U_s_doti(ix,iy,iz,V2index(iy,ix,iz)));
            U_f_V2(iy,ix,iz) = squeeze(U_f_doti(ix,iy,iz,V2index(iy,ix,iz)));
        end
    end
end
U_s_V1 = flip(flip(U_s_V1,1),2);U_f_V1 = flip(flip(U_f_V1,1),2);
U_s_V2 = flip(flip(U_s_V2,1),2);U_f_V2 = flip(flip(U_f_V2,1),2);

% montagestack(abs(U_s_V1).*mask);caxis([0 1])
% montagestack(abs(U_f_V1).*mask);caxis([0 1])
% montagestack(abs(U_s_V2).*mask);caxis([0 1])
% montagestack(abs(U_f_V2).*mask);caxis([0 1])

Atot = A1x+A2x;
A1norm = A1x./Atot.*mask;
A2norm = A2x./Atot.*mask;

U_s_V1w = abs(U_s_V1).*A1norm; %montagestack(U_s_V1w);caxis([0 1])
U_f_V1w = abs(U_f_V1).*A1norm; %montagestack(U_f_V1w);caxis([0 1])
U_s_V2w = abs(U_s_V2).*A2norm; %montagestack(U_s_V2w);caxis([0 1])
U_f_V2w = abs(U_f_V2).*A2norm; %montagestack(U_f_V2w);caxis([0 1])

U_s = U_s_V1w+U_s_V2w; U_f = U_f_V1w+U_f_V2w;
montagestack(U_s);caxis([0 1])
montagestack(U_f);caxis([0 1])

save('slowfasttwo.mat','U_s','U_f','U_s_V1','U_f_V1','U_s_V2','U_f_V2')

