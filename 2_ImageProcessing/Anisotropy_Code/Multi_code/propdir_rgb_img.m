clear all
dfile = dir('hockey*.mat');
dfile2 = dir('BMR.mat');
load(dfile.name)
load(dfile2.name)

Xsmooth = (2.446/1.764)*imgaussfilt3(real(Xmotion),1.5)+i*imgaussfilt3(imag(Xmotion),1.5);
Ysmooth = (2.446/1.764)*imgaussfilt3(real(Ymotion),1.5)+i*imgaussfilt3(imag(Ymotion),1.5);
Zsmooth = (2.446/1.764)*imgaussfilt3(real(Zmotion),1.5)+i*imgaussfilt3(imag(Zmotion),1.5);

Xg = repmat(1:120,[120 1 64])*3000;
Yg = repmat((1:120)',[1 120 64])*3000;
Zg = repmat(reshape(1:64,[1 1 64]),[120 120 1])*3000;

% h = 1/9*ones(3);
% for ii = 1:size(Xmotion,3)
% Xsmooth(:,:,ii) = filter2(h,real(Xmotion(:,:,ii)))+i*filter2(h,imag(Xmotion(:,:,ii)));
% Ysmooth(:,:,ii) = filter2(h,real(Ymotion(:,:,ii)))+i*filter2(h,imag(Ymotion(:,:,ii)));
% Zsmooth(:,:,ii) = filter2(h,real(Zmotion(:,:,ii)))+i*filter2(h,imag(Zmotion(:,:,ii)));
% end

slbottom = 1;
sltop = 64;

Xmotion = Xsmooth(:,:,slbottom:sltop);
Ymotion = Ysmooth(:,:,slbottom:sltop);
Zmotion = Zsmooth(:,:,slbottom:sltop);
mask = mask(:,:,slbottom:sltop);
Xmotion((Xmotion==0)) = NaN;
Ymotion((Ymotion==0)) = NaN;
Zmotion((Zmotion==0)) = NaN;
[Brfo_all,Bifo_all,X,Y,Z, coef,BulkMotion]=calc_3Ddisp_BulkMotionRemoved(Xmotion,Ymotion,Zmotion);

[CURLXr, CURLYr, CURLZr, ~] = curl(Xg, Yg, Zg, real(Xmotion),real(Ymotion),real(Zmotion));
[CURLXi, CURLYi, CURLZi, ~] = curl(Xg, Yg, Zg, imag(Xmotion),imag(Ymotion),imag(Zmotion));

Brfo_all = cat(4,CURLXr,CURLYr,CURLZr);
Bifo_all = cat(4,CURLXi,CURLYi,CURLZi);


i=sqrt(-1);         % imaginary unit
slice = 25;         % SLICE TO VISUALIZE
slice = slice-(slbottom-1);
usecomplexmagnitude = 1; % if 0, fit G'; if 1 fit G*
FAthresh = 0.7;     % threshold for anisotropic matter only
polthresh = 0.5;    % threshold for polarization vector orientation
ampthresh0 = 0.5;   % threshold fo amplitude = scale * median amplitude
nampthresh0 = 0.0;   % threshold fo prop vec amplitude = scale * median amplitude
kernelsize = 7;     % use MRE results with kernel size of kernelsize
nscal = 1;          % scale color map of prop dir
% U = squeeze(Brfo_all(:,:,:,1) + i *Bifo_all(:,:,:,1));
% V = squeeze(Brfo_all(:,:,:,2) + i *Bifo_all(:,:,:,2));
% W = squeeze(Brfo_all(:,:,:,3) + i *Bifo_all(:,:,:,3));
%% Analyze one component of vector displacement field at a time
for ndir = 1:3,
% CONSTRUCT DISPLACEMENT COMPONENT FIELD FROM REAL AND IMAG PARTS
% IN SPECIFIC MOTION ENCODING GRADIENT DIRECTION
U = squeeze(Brfo_all(:,:,:,ndir) + i *Bifo_all(:,:,:,ndir));
AMP(ndir) = nanstd(U(:));
%
% Get rid of NaNs for processing
I = find(isnan(U));
U(I)=zeros(size(I));
% Get Cartesian coordinates
[NUMROW,NUMCOL,NUMSLICE]=size(U);
[x,y,z]=meshgrid(1:NUMCOL,1:NUMROW,1:NUMSLICE);
% ESTIMATE PROPAGATION direction vectors
[NNX,NNY,NNZ]=dirfilt_3D_vecs_orig(U,mask);
% axis off
NNXX(:,:,:,ndir) = NNX;
NNYY(:,:,:,ndir) = NNY;
NNZZ(:,:,:,ndir) = NNZ;
UU(:,:,:,ndir) = U.*mask;
end;
% replace zeros with NaNs
UU(find(~abs(UU)))=NaN*ones(size(find(~abs(UU))));
% initialize PROPAGATION direction array for cross product calculation
NNU = zeros(NUMROW,NUMCOL,NUMSLICE,3);
% propdir will bedifferent for each MEG (u,v,w)
% re-order PROPAGATION direction for cross product calculation
NNU(:,:,:,1) = sum(NNXX,4)/3;
NNU(:,:,:,2) = sum(NNYY,4)/3;
NNU(:,:,:,3) = sum(NNZZ,4)/3;
for ii = 1:NUMROW
for jj = 1:NUMCOL
NNUnorm(ii,jj) = norm(squeeze(NNU(ii,jj,slice,:)));
end
end
nscal = 2.5/max(NNUnorm(:));
% nscal = 0.7;
% SET UP FOR PLOTTING
% for color maps
% Propagation direction vectors
% FOR RGB COLOR MAPS OF PROPAGATION DIRECTION
% NNU = permute(NNU,[3 1 2 4]);
imn(:,:,1) = nscal*NNU(:,:,slice,1);
imn(:,:,2) = nscal*NNU(:,:,slice,2);
imn(:,:,3) = nscal*NNU(:,:,slice,3);

for slice = 30:3:33
nscal = 0.75;
imn(:,:,1) = nscal*V1x(:,:,slice,1);
imn(:,:,2) = nscal*V1x(:,:,slice,2);
imn(:,:,3) = nscal*V1x(:,:,slice,3);
figure;image(abs(imn)),axis equal,axis off
% print(gcf,'-dpng')
nscal = 0.75;
imn(:,:,1) = nscal*V2x(:,:,slice,1);
imn(:,:,2) = nscal*V2x(:,:,slice,2);
imn(:,:,3) = nscal*V2x(:,:,slice,3);
 figure;image(abs(imn)),axis equal,axis off
% print(gcf,'-dpng')
% figure;imagesc(mask(:,:,slice));axis equal,axis off;colormap gray;
% print(gcf,'-dpng')
end

for slice = 18:3:45
nscal = 1;
imn(:,:,1) = nscal*N_norm(:,:,slice,1);
imn(:,:,2) = nscal*N_norm(:,:,slice,2);
imn(:,:,3) = nscal*N_norm(:,:,slice,3);
figure;image(imn),axis equal,axis off
end

