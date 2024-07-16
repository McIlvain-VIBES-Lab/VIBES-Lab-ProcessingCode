function [Ux,Uy,Uz,Theta0,Phi0] = buckyanalyzer(Bifo_all,Brfo_all,mask,cartmat)
U3 = Brfo_all + 1i.*Bifo_all;
Uo = zeros(size(U3));
xy = length(U3(:,1,1,1));
z = length(U3(1,1,:,1));

for dir = 1:3
tic
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
    KX=KX-cenx; KY=KY-ceny; KZ=KZ-cenz;        
    % magnitude of wavenumber if any radial filtering
    KR = sqrt(KX.^2+KY.^2+KZ.^2);        

    %load buckyinfo.mat allnvec;
    allnvec = cartmat;
    Theta0 = zeros(size(allnvec(:,1)));
    Phi0 = zeros(size(allnvec(:,1)));
    for nn=1:length(allnvec);        
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
        Gf = G3.*cn2;
        %   Uf(:,:,:,nn)=Gf;
        % Uf(:,:,:,nn) = (ifftn(ifftshift(Gf)));
        Ui(:,:,:,nn) = (ifftn(ifftshift(Gf)));
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if dir == 1
        Ux = Ui;
        'Ux found'
    elseif dir == 2;
        Uy = Ui;
        'Uy found'
    else 
        Uz = Ui;
        'Uz found'
    end
toc
end
