function [Ui]=dirfilt_3D_vecs_alldirs_v2(Brfo_all,Bifo_all,mask,A)
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
toc
sprintf('U found')
tic
load 'allnvec_matrix.mat' allnvec_matrix 
%load in matrix for crossing with A, the fiber direction from DTI
allnvec_matrix(abs(allnvec_matrix)<=.005) = 0; %matrix isn't zeroed out for some reason so just sets numbers to 0 that should be 0
m_s = zeros(80,80,48,3,92);
m_f = zeros(size(m_s));
%Individual cross products for each of the 92 directions for m_s and m_f
% m_s = NxA
% m_f = m_sxA
% for nn = 1:92
%     m_s(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),A,4); 
%     m_f(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),squeeze(m_s(:,:,:,:,nn)),4);
% end
 for nn = 1:92
     m_s(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),A,4); 
     m_f(:,:,:,:,nn) = cross(squeeze(allnvec_matrix(:,:,:,:,nn)),squeeze(m_s(:,:,:,:,nn)),4);
 end
toc
sprintf('Waves found')
tic

%% take dot products and multiply by m_s and m_f to create U vectors
U_s = zeros(size(m_s));
U_f = zeros(size(U_s)); 
U_s_dot = zeros(80,80,48,92);
U_f_dot = zeros(size(U_s_dot));
for nn = 1:length(allnvec(:,1))
    %individual dot products for U_s and U_f
    U_s_dot(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_s(:,:,:,:,nn)),4)); 
    U_f_dot(:,:,:,nn) = squeeze(dot(Ui(:,:,:,:,nn),squeeze(m_f(:,:,:,:,nn)),4));
    %taking the individual dot products and then using it as a scaling
    %measure for the individual directionality vector
    U_s(:,:,:,1,nn) = squeeze(m_s(:,:,:,1,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_s(:,:,:,:,nn)),4);
    U_s(:,:,:,2,nn) = squeeze(m_s(:,:,:,2,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_s(:,:,:,:,nn)),4);
    U_s(:,:,:,3,nn) = squeeze(m_s(:,:,:,3,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_s(:,:,:,:,nn)),4);
    U_f(:,:,:,1,nn) = squeeze(m_f(:,:,:,1,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_f(:,:,:,:,nn)),4);
    U_f(:,:,:,2,nn) = squeeze(m_f(:,:,:,2,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_f(:,:,:,:,nn)),4);
    U_f(:,:,:,3,nn) = squeeze(m_f(:,:,:,3,nn)).*dot(Ui(:,:,:,:,nn),squeeze(m_f(:,:,:,:,nn)),4);
end
toc
sprintf('individual slow and fast wave fields found')

%% Summation of the dot products and creation of N_bar vector for finding theta
%Usd_bar = sum(U_s_dot,4); Ufd_bar = sum(U_f_dot,4);
Usd_bar = sum(Us_dot,4); % Summing the dot products of the U_s (making 80x80x48)
Ufd_bar = sum(Uf_dot,4); % Summing the dot products of the U_f (making 80x80x48)

Usd_norm = abs(Usd_bar)./sqrt((Usd_bar.^2)+(Ufd_bar.^2));
Ufd_norm = abs(Ufd_bar)./sqrt((Usd_bar.^2)+(Ufd_bar.^2));

% Ucomb = sqrt(U_s_dot.^2+U_f_dot.^2); % Magnitude of the summation of the dot products
% Ucomb(Ucomb==0) = NaN; 
% Utot = cat(5,Ucomb,Ucomb,Ucomb); % Creating a 80x80x48x92x3 matrix
% Utot = permute(Utot,[1 2 3 5 4]); % changing the matrix to same size as N
% U_all = sum(Utot,5); % summing across the 92 directions

Ui_tot = sqrt(sum(Ui.^2,5));
for nn = 1:92
    Nweight(:,:,:,:,nn) = (abs(Ui(:,:,:,:,nn)).*allnvec_matrix(:,:,:,:,nn))./Ui_tot;
end
Nfinal = sum(Nweight,5);

% Ntot = zeros(size(Ui));
% for nn = 1:92
%     Ntot(:,:,:,:,nn) = Utot(:,:,:,:,nn).*Ui(:,:,:,:,nn); 
%     % Imposing the combined dot product matrix upon the original
%     % directional matrix 
% end
% U_s_all = sum(U_s,5); U_f_all = sum(U_f,5); m_s_all = sum(m_s,5); m_f_all = sum(m_f,5); 
% %U_s_norm = U_s_all./U_all; U_f_norm = U_f_all./U_all; 
% Ntot_bar = sum(Ntot,5); % Summation of the dot products across all 92x3 directions
% Usd_norm = Usd_bar./sum(Ucomb,4); % Normalized dot product for slow waves
% Ufd_norm = Ufd_bar./sum(Ucomb,4); % Normalized dot product for fast waves
% Nfinal = Ntot_bar./U_all; % Normalization of the total directionality matrix
sprintf('Solving for Theta...')

%% Theta definition in degrees
theta = abs(acosd(dot(Nfinal,A,4))); % solving for Theta (cos(theta) = N*A)
% creating the range of of the angles to 0-360
tmp = ones(size(theta(abs(theta)>360)))*360; 
theta(abs(theta)>360) = theta(abs(theta)>360)-tmp;
theta(isnan(theta)) = 0;

%% Save 

sprintf('saving...')

% save('SlowFastWaves.mat','U_s_all','U_f_all','m_s_all','m_f_all','m_s','m_f','U_s','U_f')
% save('Waveprop_BMR.mat','Ntot_bar','Usd_bar','Ufd_bar','Usd_norm','Ufd_norm','Nfinal','theta','U_all');
% save('SlowFastWaves.mat','U_s_all','U_f_all','m_s_all','m_f_all','m_s','m_f','U_s','U_f')
save('Waveprop_BMR.mat','Usd_bar','Ufd_bar','Usd_norm','Ufd_norm','Nfinal','theta');
end
