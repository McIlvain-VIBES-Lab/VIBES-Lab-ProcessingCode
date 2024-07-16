function [Uf]=dirfilt_3D_vecs_alldirs(Xmotion,Ymotion,Zmotion,mask)
% show amplitude weighted direction of wave propagation by directional
% filtering
%  % uses 92 directions based on bucky ball vertices and mid faces
% written by R. Okamoto, Washington University 2014
% inputs:
% U3 - (NxMxP) field of complex coefficients (Re,Im) of fundamental component
% mask - mask array (NxM)

% outputs
% UUX, UUY, UUZ: amplitude-weighted vectors of propagation direction

U3 = cat(4,Xmotion(:,:,:,1),Ymotion(:,:,:,2),Zmotion(:,:,:,3));

% G3 = spatial FFT of complex coefficient field, needs to be fftshifted to center it
G3 = fftshift(fftn(U3.*mask));

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
for nn=1:length(allnvec);
    
    
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
    Uf(:,:,:,nn) = (ifftn(ifftshift(Gf)));
    %%%%%%%%%%%%%%%%%%%%%%%%%
end


