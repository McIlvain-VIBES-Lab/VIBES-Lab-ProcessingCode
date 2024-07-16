function [Ux,Uy,Uz,Theta0,Phi0,U3filt] = buckyanalyzer_CLJfilter(Bifo_all,Brfo_all,mask,cartmat)
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
G3 = zeros(80,80,48,3);
U3filt = zeros(size(G3));
for dir = 1:3
    G3(:,:,:,dir) = fftshift(fftn(U3(:,:,:,dir).*mask));
    G3(:,:,:,dir) = squeeze(G3(:,:,:,dir)).*filter;
    U3filt(:,:,:,dir) = (ifftn(ifftshift(squeeze(G3(:,:,:,dir)))));
end
toc
tic
disp('G3 found')
   

% load buckyinfo.mat allnvec;
allnvec = cartmat;
Theta0 = zeros(size(allnvec(:,1)));
Phi0 = zeros(size(allnvec(:,1)));
disp('setup complete')
cn2_nn = zeros(80,80,48,size(allnvec,1));
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
    Ui = zeros(80,80,48);
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
    
%     if dir == 1
%         Ux = Ui;
%         'Ux found'
%     elseif dir == 2;
%         Uy = Ui;
%         'Uy found'
%     else 
%         Uz = Ui;
%         'Uz found'
%     end
% toc
% end
