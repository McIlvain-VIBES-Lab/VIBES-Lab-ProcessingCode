function [Ui,Ns_norm,Nf_norm,U_si,U_fi]=dirfilt_3D_vecs_alldirs_v4(Brfo_all,Bifo_all,mask,A)
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
U3 = Brfo_all + 1i.*Bifo_all;
Uo = zeros(size(U3));
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
    load buckyinfo.mat allnvec;
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
    end      
end 
toc

%% Cross products to create m_s & m_f
% N = Uf; %defines the base directional vectors
sprintf('U found')
tic
load 'allnvec_matrix.mat' allnvec_matrix 
%load in matrix for crossing with A, the fiber direction from DTI
allnvec_matrix(abs(allnvec_matrix)<=.005) = 0; %matrix isn't zeroed out for some reason so just sets numbers to 0 that should be 0

m_si = zeros(80,80,48,3,92);
m_fi = zeros(size(m_si));

m_sx = zeros(size(m_si(:,:,:,1))); m_sy = zeros(size(m_si(:,:,:,1))); m_sz = zeros(size(m_si(:,:,:,1)));
m_si_unit = zeros(size(m_si));

m_fx = zeros(size(m_si(:,:,:,1))); m_fy = zeros(size(m_si(:,:,:,1))); m_fz = zeros(size(m_si(:,:,:,1)));
m_fi_unit = zeros(size(m_si));

%Individual cross products for each of the 92 directions for m_s and m_f
% m_s = NxA
% m_f = m_sxA

 for nn = 1:92
     m_si(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),A,4); 
     m_sx(:,:,:,nn) = squeeze(m_si(:,:,:,1,nn)); m_sy(:,:,:,nn) = squeeze(m_si(:,:,:,2,nn)); m_sz(:,:,:,nn) = squeeze(m_si(:,:,:,3,nn)); 
     m_si_unit(:,:,:,1,nn) = m_sx(:,:,:,nn)./(sqrt(m_sx(:,:,:,nn).^2+m_sy(:,:,:,nn).^2+m_sz(:,:,:,nn).^2));
     m_si_unit(:,:,:,2,nn) = m_sy(:,:,:,nn)./(sqrt(m_sx(:,:,:,nn).^2+m_sy(:,:,:,nn).^2+m_sz(:,:,:,nn).^2));
     m_si_unit(:,:,:,3,nn) = m_sz(:,:,:,nn)./(sqrt(m_sx(:,:,:,nn).^2+m_sy(:,:,:,nn).^2+m_sz(:,:,:,nn).^2));
     
     m_fi(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),squeeze(m_si_unit(:,:,:,:,nn)),4);
     m_fx(:,:,:,nn) = squeeze(m_fi(:,:,:,1,nn)); m_fy(:,:,:,nn) = squeeze(m_fi(:,:,:,2,nn)); m_fz(:,:,:,nn) = squeeze(m_fi(:,:,:,3,nn)); 
     m_fi_unit(:,:,:,1,nn) = m_fx(:,:,:,nn)./(sqrt(m_fx(:,:,:,nn).^2+m_fy(:,:,:,nn).^2+m_fz(:,:,:,nn).^2));
     m_fi_unit(:,:,:,2,nn) = m_fy(:,:,:,nn)./(sqrt(m_fx(:,:,:,nn).^2+m_fy(:,:,:,nn).^2+m_fz(:,:,:,nn).^2));
     m_fi_unit(:,:,:,3,nn) = m_fz(:,:,:,nn)./(sqrt(m_fx(:,:,:,nn).^2+m_fy(:,:,:,nn).^2+m_fz(:,:,:,nn).^2));
 end
toc
sprintf('Waves found')
tic

%% take dot products and multiply by m_s and m_f to create U vectors
U_si = zeros(size(m_si_unit));
U_fi = zeros(size(U_si)); 
U_s_doti = zeros(size(squeeze(m_si_unit(:,:,:,1,:))));
U_f_doti = zeros(size(U_s_doti));
for nn = 1:length(allnvec(:,1))
    %individual dot products for U_s and U_f
    U_s_doti(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_si_unit(:,:,:,:,nn)),4)); 
    U_f_doti(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_fi_unit(:,:,:,:,nn)),4));
    %taking the individual dot products and then using it as a scaling
    %measure for the individual directionality vector
    for dir = 1:3
        U_si(:,:,:,dir,nn) = squeeze(m_si_unit(:,:,:,dir,nn)).*U_s_doti(:,:,:,nn);
        U_fi(:,:,:,dir,nn) = squeeze(m_fi_unit(:,:,:,dir,nn)).*U_f_doti(:,:,:,nn);
    end
end
toc
sprintf('individual slow and fast wave fields found')

%% Summation of the dot products and creation of N_bar vector for finding theta
U_s = sum(abs(U_si),5); U_f = sum(abs(U_fi),5);
Us_dot = sum(abs(U_s_doti),4); Uf_dot = sum(abs(U_f_doti),4);
U_dot_all = Us_dot+Uf_dot;

U_sfi = U_si+U_fi; 
U_sf = U_s+U_f;

for nn = 1:92
    Ntot(:,:,:,:,nn) = U_sfi(:,:,:,:,nn).*allnvec_matrix(:,:,:,:,nn); 
    % Imposing the combined dot product matrix upon the original
    % directional matrix 
end

Nweighted = sum(abs(Ntot),5).*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
Nfinal = Nweighted./U_sf.*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
Nfinx = Nfinal(:,:,:,1); Nfiny = Nfinal(:,:,:,2); Nfinz = Nfinal(:,:,:,3);
Nfinal_dir(:,:,:,1) = Nfinx./(sqrt(Nfinx.^2+Nfiny.^2+Nfinz.^2));
Nfinal_dir(:,:,:,2) = Nfiny./(sqrt(Nfinx.^2+Nfiny.^2+Nfinz.^2));
Nfinal_dir(:,:,:,3) = Nfinz./(sqrt(Nfinx.^2+Nfiny.^2+Nfinz.^2));

Us_norm_i = abs(U_si)./abs(U_sfi); Uf_norm_i = abs(U_fi)./abs(U_sfi);
Ns_norm_i = Us_norm_i.*allnvec_matrix; Nf_norm_i = Uf_norm_i.*allnvec_matrix;

Us_norm_i(isnan(Us_norm_i)) = 0; Uf_norm_i(isnan(Uf_norm_i)) = 0;
Ns_norm_i(isnan(Ns_norm_i)) = 0; Nf_norm_i(isnan(Nf_norm_i)) = 0;

m_s = mean(m_si_unit,5); m_f = mean(m_fi_unit,5);
Us_norm = mean(abs(Us_norm_i),5).*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
Uf_norm = mean(abs(Uf_norm_i),5).*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
Ns_norm = mean(abs(Ns_norm_i),5).*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
Nf_norm = mean(abs(Nf_norm_i),5).*repmat(flip(flip(permute(mask,[2 1 3]),2),1),[1 1 1 3]);
%Ns_norm(Ns_norm==0) = NaN; %Nf_norm(Nf_norm==0) = NaN;

Ns_vox = Ns_norm>Nf_norm; Nf_vox = Nf_norm>Ns_norm;

S_vox_U = abs(sum(Us_norm,4)>sum(Uf_norm,4)); F_vox_U = abs(sum(Uf_norm,4)>sum(Us_norm,4)); 
F_vox = abs(sum(Nf_norm,4)>sum(Ns_norm,4)); S_vox = sum(Ns_norm,4)>sum(Nf_norm,4);

%% Save 

status = mkdir('Iterative_Outputs');

cd Iterative_Outputs
sprintf('saving...')

%save('SlowFastWaves.mat','U_si','U_fi','m_s','m_f','U_s','U_f')
%save('Waveprop_BMR.mat','Us_dot','Uf_dot','Nfinal','Ns_norm','Nf_norm','U_dot_all','Us_norm','Uf_norm','Nweighted','Nfinal_dir');
%save('Outputs.mat','Ns_vox','Nf_vox','F_vox','S_vox','S_vox_U','F_vox_U')
cd ..
end
