function [Gp, Gdp,C, nre, R, dnre,smoothedD,smoothedL]=fitsmooth_uvw_rjofun_spiral_R2_parfor(Brfo_all,Bifo_all,mask,voxsize, PC2mm, fms, rho, h,Delta,nanfract,Nxy, usecurl,curldata,smoothd,smoothL,fname) 
% function to fit 3-D MRE data using TLS
% inputs:
%Brfo_all,
%Bifo_all,
%mask,
%voxsize, 
%PC2mm, 
%fms, 
%rho,
%h
%Delta,
%nanfract
%Nxy
%usecurl
%smoothd
%smoothL
%fname
% curldata
% written by R. Okamoto  2013-08-05
% RJO revised 2015-06-06 to use polynomial curl
% RJO revised 2015-09-08 to write out data directly to file
% RJO revised 2015-09-25 to write out fms, h, PC2mm
% RJO revised 2016-10-25 to run loops in parallel
% RJO revised 2016-11-03 to put smoothedL and smoothD and smoothL in the output file
% 
%% 
% fms=100;
% h =2e-3;      % 2 mm voxel
if exist('smoothd','var'),
    if ~isempty(smoothd);
        sd = smoothd(1);    % standard deviation for Gaussian smoothing in pix
        convsize=smoothd(2); % size of convolution kernel for Gaussian smoothing, corners have 3.4% weight of center
        
    end
end
if exist('smoothL','var')% then smooth laplacian
    if ~isempty(smoothL);
        sdL = smoothL(1);
        convsizeL=smoothL(2);
    end
end
%    umask = 0; % mask before smoothing?
%Delta = 1;
%nanfract = 0.75;

%% get size of image
[npy, npx, nslice,nMEG] = size(Brfo_all);
% create a grid of the right dimensions
[x,y,z] = meshgrid(1:1:npx,1:1:npy,-Nxy:1:(nslice+Nxy-1));
% scale values to be mm instead of voxels
x = x*voxsize(1); y = y*voxsize(2); z = z*voxsize(3);

%% pad displacement and mask with NaN
Brfo_pad(:,:,Nxy+(1:nslice),:)=Brfo_all;
Brfo_pad(:,:,[1:Nxy Nxy+nslice+(1:Nxy)],:)=NaN;
Bifo_pad(:,:,Nxy+(1:nslice),:)=Bifo_all;
Bifo_pad(:,:,[1:Nxy Nxy+nslice+(1:Nxy)],:)=NaN;
mask_pad(:,:,Nxy+(1:nslice))=mask;
mask_pad(:,:,[1:Nxy Nxy+nslice+(1:Nxy)])=NaN;
%% define displacement components,mask and pad with NaN 
u=PC2mm*(Brfo_pad(:,:,:,1)+1i*Bifo_pad(:,:,:,1)).*mask_pad;
v=PC2mm*(Brfo_pad(:,:,:,2)+1i*Bifo_pad(:,:,:,2)).*mask_pad;
w=PC2mm*(Brfo_pad(:,:,:,3)+1i*Bifo_pad(:,:,:,3)).*mask_pad;
%% 
if usecurl, % take curl
    %[curlx,curly,curlz,cav]=curl(x,y,z,u,v,w);
    curldata_pad(:,:,Nxy+(1:nslice),:)=curldata;
    curldata_pad(:,:,[1:Nxy Nxy+nslice+(1:Nxy)],:)=NaN+NaN*1i;
    curldata_pad=real(curldata_pad).*repmat(mask_pad,[1 1 1 3])+1i*imag(curldata_pad).*repmat(mask_pad,[1 1 1 3]);
    curlx=curldata_pad(:,:,:,1);curly=curldata_pad(:,:,:,2);curlz=curldata_pad(:,:,:,3);
   

    if exist('sd','var'), %only defined if smoothing displacements
        padsz=(convsize-1)/2;
        disp(['Smoothing curl w/ kernel: ',int2str(convsize),' and SD: ',num2str(sd)])
%         curlx=  padarray(curldata(:,:,:,1),[padsz padsz padsz],NaN+1i*NaN);
%         curly=  padarray(curldata(:,:,:,2),[padsz padsz padsz],NaN+1i*NaN);
%         curlz=  padarray(curldata(:,:,:,3),[padsz padsz padsz],NaN+1i*NaN);
        
        curlx = smooth3nan(curlx,'gaussian',convsize,sd);
        curly = smooth3nan(curly,'gaussian',convsize,sd);
        curlz = smooth3nan(curlz,'gaussian',convsize,sd);
        
%         curlx=curlx((1+padsz):(end-padsz),(1+padsz):(end-padsz),(1+padsz):(end-padsz));
%         curly=curly((1+padsz):(end-padsz),(1+padsz):(end-padsz),(1+padsz):(end-padsz));
%         curlz=curlz((1+padsz):(end-padsz),(1+padsz):(end-padsz),(1+padsz):(end-padsz));

    end
    smoothedD(:,:,:,1)=curlx;
    smoothedD(:,:,:,2)=curly;
    smoothedD(:,:,:,3)=curlz;

else
    if exist('smoothd'); % smooth displacement field
        u=smooth3nan(u,'gaussian',convsize,sd);
        v=smooth3nan(v,'gaussian',convsize,sd);
        w=smooth3nan(w,'gaussian',convsize,sd);
    end
    smoothedD(:,:,:,1)=u;
    smoothedD(:,:,:,2)=v;
    smoothedD(:,:,:,3)=w;

end


% take derivatives
if usecurl,
    Lcx = 6*del2(curlx,h);
    Lcy = 6*del2(curly,h);
    Lcz = 6*del2(curlz,h);
    if exist('smoothL'),
        if ~isempty(smoothL);
            disp(['Smoothing Laplacian of curl ...'])
            Lcx = smooth3nan(Lcx,sdL,convsizeL);
            Lcy = smooth3nan(Lcy,sdL,convsizeL);
            Lcz = smooth3nan(Lcz,sdL,convsizeL);     
        end
    end
    smoothedL(:,:,:,1)=Lcx;
    smoothedL(:,:,:,2)=Lcy;
    smoothedL(:,:,:,3)=Lcz;
    
    % mask curl and laplacian
%         curlx=curlx.*mask;
%         curly=curly.*mask;
%         curlz=curlz.*mask;
%         
%         Lcx=Lcx.*mask;
%         Lcy=Lcy.*mask;
%         Lcz=Lcz.*mask;
        

else % use displacements
    Lu = 6*del2(u,h);
    Lv = 6*del2(v,h);
    Lw = 6*del2(w,h);
    if exist('smoothL'),
        if ~isempty(smoothL),
        disp(['smoothing Laplacian of displacement ...']);
        Lu = smooth3nan(Lu,sdL,convsizeL);
        Lv = smooth3nan(Lv,sdL,convsizeL);
        Lw = smooth3nan(Lw,sdL,convsizeL);
        end
    end
    
    smoothedL(:,:,:,1)=Lu;
    smoothedL(:,:,:,2)=Lv;
    smoothedL(:,:,:,3)=Lw;
%if exist (smoothL)
    % mask Laplacian
%     Lu=Lu.*mask;
%     Lv=Lv.*mask;
%     Lw=Lw.*mask;
end; %if usecurl

% loop through  the slices
allslices  = 1:nslice;
disp(['Using slices ',int2str(allslices(1)),' to ',int2str(allslices(end))])

poolobj = parpool; % 
if usecurl,
    parfor k=1:(length(allslices));
    slices=(k+(0:(2*Nxy)));

    disp(['Fitting slices ',int2str(slices(1)-Nxy),' to ',int2str(slices(end)-Nxy), ' with Delta = ',num2str(Delta)])
 %   [Gp(:,:,k), Gdp(:,:,k), C(:,:,:,k), nre(:,:,k), R(:,:,k), dnre(:,:,k)] = matLSfit3DMS_v3d(uu,vv,ww,del2uu,del2vv,del2ww,rho,fms,Nxy,maskk,Delta,nanfract);
  [Gp(:,:,k), Gdp(:,:,k), C(:,:,:,k), nre(:,:,k), R(:,:,k), dnre(:,:,k)] =matLSfit3DMS_v3d(curlx(:,:,slices),curly(:,:,slices),curlz(:,:,slices),Lcx(:,:,slices),Lcy(:,:,slices),Lcz(:,:,slices),rho,fms,Nxy,mask_pad(:,:,slices),Delta,nanfract);
    end  %for k (par)
else
   parfor k=1:(length(allslices));
    slices=(k+(0:(2*Nxy)));
   [Gp(:,:,k), Gdp(:,:,k), C(:,:,:,k), nre(:,:,k), R(:,:,k), dnre(:,:,k)] =matLSfit3DMS_v3d(u(:,:,slices),v(:,:,slices),w(:,:,slices),Lu(:,:,slices),Lv(:,:,slices),Lw(:,:,slices),rho,fms,Nxy,mask_pad(:,:,slices),Delta,nanfract);
   end %parfor
end % if
delete(poolobj)   % close pool
if exist('fname','var'),
    save(fname,'Gp','Gdp','nre','smoothedD','smoothedL','Delta','nanfract','Nxy','usecurl','mask','fms','PC2mm','h','smoothd','smoothL');
end
