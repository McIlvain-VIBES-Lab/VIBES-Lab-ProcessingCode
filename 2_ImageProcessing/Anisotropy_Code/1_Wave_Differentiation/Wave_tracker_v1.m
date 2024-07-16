function Wave_tracker_v1()
%% Pick starting point for seed
% Identify the edge of the mask
% es = 0;
% for ii = 1:size(mask,3)
%     if any(any(mask(:,ii,:)))
%         if es == 0
%             es = ii;
%         else
%         end
%     else
%     end
% end
% 
% % Find centroid of this side (starting point)
% cent = regionprops(squeeze(mask(:,es,:)));
% cent1 = round(cent.Centroid,0);
% xx = cent1(2);
% zz = cent1(1);
% yy = es;

%% Find Primary Direction for seed

%Filter
x = length(Brfo_all(1,:,1,1));
y = length(Brfo_all(:,1,1,1));
z = length(Brfo_all(1,1,:,1));

Xsmooth = imgaussfilt3(Brfo_all(:,:,:,1),1.5)+1i*imgaussfilt3(Bifo_all(:,:,:,1),1.5);
Ysmooth = imgaussfilt3(Brfo_all(:,:,:,2),1.5)+1i*imgaussfilt3(Bifo_all(:,:,:,2),1.5);
Zsmooth = imgaussfilt3(Brfo_all(:,:,:,3),1.5)+1i*imgaussfilt3(Bifo_all(:,:,:,3),1.5);

Xg = repmat(1:x,[y 1 z])*3000;
Yg = repmat((1:y)',[1 x z])*3000;
Zg = repmat(reshape(1:z,[1 1 z]),[y x 1])*3000;

slbottom = 1;
sltop = z;

Xmotion = Xsmooth(:,:,slbottom:sltop);
Ymotion = Ysmooth(:,:,slbottom:sltop);
Zmotion = Zsmooth(:,:,slbottom:sltop);
mask = mask(:,:,slbottom:sltop);
Xmotion((Xmotion==0)) = NaN;
Ymotion((Ymotion==0)) = NaN;
Zmotion((Zmotion==0)) = NaN;

[CURLXr, CURLYr, CURLZr, ~] = curl(Xg, Yg, Zg, real(Xmotion),real(Ymotion),real(Zmotion));
[CURLXi, CURLYi, CURLZi, ~] = curl(Xg, Yg, Zg, imag(Xmotion),imag(Ymotion),imag(Zmotion));

Brfo_all = cat(4,CURLXr,CURLYr,CURLZr);
Bifo_all = cat(4,CURLXi,CURLYi,CURLZi);

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
filter(((xy/2)+1),((xy/2)+1),(ceil((z/2)+1))) = 0;

tic
G3 = zeros(y,x,z,3);
U3filt = zeros(size(G3));
for dir = 1:3
    G3(:,:,:,dir) = fftshift(fftn(U3(:,:,:,dir).*mask));
    G3(:,:,:,dir) = squeeze(G3(:,:,:,dir)).*filter;
    U3filt(:,:,:,dir) = (ifftn(ifftshift(squeeze(G3(:,:,:,dir)))));
end
toc
tic
disp('G3 found')
   

load bucky200.mat allnvec;
%load bucky150.mat allnvec
%allnvec = X;
Theta0 = zeros(size(allnvec(:,1)));
Phi0 = zeros(size(allnvec(:,1)));
disp('setup complete')
cn2_nn = zeros(y,x,z,size(allnvec,1));
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
    Ui = zeros(y,x,z);
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

tic
Ut = cat(5,Ux,Uy,Uz);
UTa = sqrt(sum(abs(Ut).^2,5));
U3a = sqrt(sum(abs(U3filt).^2,4));
UT_norm = UTa./(repmat(U3a,[1 1 1 size(UTa,4)]));
UT_mask = UT_norm(:,:,:,:).*repmat(mask(:,:,:),[1 1 1 size(UTa,4)]);
toc

%% Seed Point Identification

[SeedPoints] = SeedPointFinder(UT_mask,mask,allnvec);


X1_all = repmat(zeros(size(mask)),[1 1 1 length(SeedPoints)]);
XTP_all = repmat(zeros(size(mask)),[1 1 1 length(SeedPoints)]);
XAmp_all = repmat(zeros(size(mask)),[1 1 1 length(SeedPoints)]);
Xnvec_all =repmat(zeros(size(mask)),[1 1 1 3 length(SeedPoints)]);

for ss = 1:length(SeedPoints)
    ss
    xx = SeedPoints(ss,1);
    yy = SeedPoints(ss,2);
    zz = SeedPoints(ss,3);
%%  Peak Identification for seed point
% points around peak = 7
% ptx = 42;
filtn = length(allnvec);

% peak 1
tmp = squeeze(UT_mask(xx,yy,zz,:)); % (squeeze(UT_mask(37,26,28,:)));
ix = find(tmp==max(tmp));

for ii = 1:length(allnvec)
    dotv(ii,1) = dot(allnvec(ix,:),allnvec(ii,:));
end

dtheta = acosd(dotv);

[~,I] = sort(dtheta);
dfx = cat(2,I,dtheta(I),tmp(I));
ptx = find(dfx(:,2)>=50,1);

tmpx = tmp(I);
allnvecx = allnvec(I,:);

wv = allnvecx(1:ptx,:).*repmat(tmpx(1:ptx),[1 3]);
nvec = mean(wv,1)./norm(mean(wv,1));

Amp = tmpx(1); % Amp = norm(mean(wv,1));

%% SETUP for wave finding loop
[Surf_mask] = Surface_masker(mask);
X1 = ones(size(mask))-mask+Surf_mask; % X2 = ones(size(mask))-mask-Surf_mask;
% Y1 = ones(size(mask))-mask; Y2 = ones(size(mask))-mask;
% Z1 = ones(size(mask))-mask; Z2 = ones(size(mask))-mask;

X_TP = zeros(size(X1));
X_Amp = zeros(size(X1));
X_nvec = zeros(size(repmat(X1,[1 1 1 3])));

X_nvec(xx,yy,zz,:) = nvec;
X_Amp(xx,yy,zz) = Amp; % Change variable name
X1(xx,yy,zz) = 2; %Begin building area of analysis (change Variable Name)
theta_prime_max = 20;

ref_list = [xx,yy,zz];

%% WAVE FINDING LOOP

while size(ref_list,1) >= 1
    xi = ref_list(1,1);
    yi = ref_list(1,2);
    zi = ref_list(1,3);
    for ii = xi-1:xi+1
        for jj = yi-1:yi+1
            for kk = zi-1:zi+1
                if kk > size(X1,3)
                    kk = size(X1,3);
                elseif kk < 1
                    kk = 1;
                else
                end
                if X1(ii,jj,kk) == 0
                    % X1(ii,jj,kk) = 2;
                    ref_list = [ref_list;[ii,jj,kk]];
                    [nvec_new, theta_prime, Amp_new, output] = angle_compare(ii, jj, kk, nvec, theta_prime_max, allnvec, UT_mask);
                    X1(ii,jj,kk) = output;
                    X_TP(ii,jj,kk) = theta_prime; % Change variable name
                    X_Amp(ii,jj,kk) = Amp_new; %Change variable name 
                    X_nvec(ii,jj,kk,:) = nvec_new; %Change Variable nam
                else
                end
            end
        end
    end
    ref_list(1,:) = [];
end

X1_all(:,:,:,ss) = X1;
XTP_all(:,:,:,ss) = X_TP;
XAmp_all(:,:,:,ss) = X_Amp;
Xnvec_all(:,:,:,:,ss) = X_nvec;

end

save('Wave_Tracking.mat','X1_all','XTP_all','XAmp_all','Xnvec_all','SeedPoints')

%% For all connecting voxels, find most similiar peak(minimial angle between two peaks )

%% Label whether wave compoenent has been found for current wave