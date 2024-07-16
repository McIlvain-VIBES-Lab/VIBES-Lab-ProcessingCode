function [U,m_s,m_f,M_S,M_F] = test(A,Xmotion,Ymotion,Zmotion,dir,mask)

U3(:,:,:,1) = Xmotion; U3(:,:,:,2) = Ymotion; U3(:,:,:,3) = Zmotion;

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
    %load buckyinfo.mat allnvec;
    %for nn=1:length(allnvec);
     load buckyinfo.mat allnvec;
    nn_max = 5;
    for nn=1:nn_max;    
        
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
        
        % vector fields from weighted amplitudes multiplied by direction
        ux = squeeze(abs(Uf(:,:,:,nn)))*no0(1);
        uy = squeeze(abs(Uf(:,:,:,nn)))*no0(2);
        uz = squeeze(abs(Uf(:,:,:,nn)))*no0(3);
%         UX(:,:,:,nn) = ux;       % x component
%         UY(:,:,:,nn) = uy;        % y component
%         UZ(:,:,:,nn) = uz;        % z component
        U = cat(4,ux,uy,uz);
        %end;
        m_s = zeros(80,80,48,3,nn_max);
        m_f = zeros(size(m_s));
%         size(U)
%         size(A)
        for i = 1:80
            for j = 1:80
                for k = 1:48
                    tmp1 = squeeze(U(i,j,k,:));
                    tmp2 = squeeze(A(i,j,k,:));
                    m_s(i,j,k,:,nn) = cross(tmp1,tmp2);
                    m_ss = squeeze(m_s(i,j,k,:,nn));
                    m_f(i,j,k,:,nn) = cross(tmp1,m_ss);
                end
            end
        end
    end
    M_S = zeros(80,80,48,3);
    M_F = zeros(size(M_S));
%     for i = 1:80
%         for j = 1:80
%             for k = 1:48
%                 for dir1 = 1:3
%                     M_S(i,j,k,dir1) = squeeze(mean(m_s(i,j,k,dir1,:),5));
%                     M_F(i,j,k,dir1) = squeeze(mean(m_f(i,j,k,dir1,:),5));
%                 end
%             end
%         end
%     end
     for i = 1:80
         for j = 1:80
             for k = 1:48
                 for dir1 = 1:3
                     M_S_max(i,j,k,dir1) = squeeze((m_s(i,j,k,dir1,:)));
                     M_F_max(i,j,k,dir1) = squeeze((m_f(i,j,k,dir1,:)));
                 end
             end
         end
     end
    % weighted average direction of propagation
    %size(mask)
    %size(UUX(:,:,:,dir))
    %size(sum(UX,4))
%     UUX(:,:,:,dir)=-sum(UX,4).*mask;  %sum over all the filter directions
%     UUY(:,:,:,dir)=-sum(UY,4).*mask;
%     UUZ(:,:,:,dir)=-sum(UZ,4).*mask;
