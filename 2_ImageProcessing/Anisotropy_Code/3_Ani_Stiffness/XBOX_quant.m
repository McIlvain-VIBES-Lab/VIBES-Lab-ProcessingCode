load adata_xbox.mat
DTI = abs(avec_img(:,:,:,1))>0;
% cd Buckley_Tbi_Brain
cd XBOXXY
load XBox_SetupXDirY150_sys2.mat
maskXY = flip(flip(permute(mask,[2 1 3]),1),2);
cd ..; cd XBOXXZ
load XBox_SetupXDirZ150_sys2.mat
maskXZ = flip(flip(permute(mask,[2 1 3]),1),2);
cd ..;

FAtractlistAP = [];
FAtractlistLR = [];
clear ThlistLR ThlistAP MutotlistLR MutotlistAP TotalLR TotalAP AmplistAP AmplistLR UperctotAP UperctotLR

%%  AP
FA = FA_img;
%load SlowFastSummation.mat
cd XBOXXY
%addpath('PostInv')
load Mu_NLI.mat; load slowfasttwo.mat; load BuckyOut.mat;

%Vmult1 = flip(flip(permute(Vmult,[2 1 3]),1),2);
Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
Theta1 = acos(dot(avec_img,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(avec_img,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));
A1 = flip(flip(permute(A1x,[2 1 3]),1),2); 
A2 = flip(flip(permute(A2x,[2 1 3]),1),2);
UsV1 = flip(flip(permute(U_s_V1,[2 1 3]),1),2); 
UfV1 = flip(flip(permute(U_f_V1,[2 1 3]),1),2);
UsV2 = flip(flip(permute(U_s_V2,[2 1 3]),1),2); 
UfV2 = flip(flip(permute(U_f_V2,[2 1 3]),1),2);

maskx = DTI.*maskXY.*maskXZ;
list = maskx(maskx==1);

%% Mu

nn = Mux.*maskx; nn(isnan(nn))  = 0;
Mulist = nn(nn>0);
figure;im(nn);caxis([0 5000])
MutotlistAP(1:length(Mulist)) = Mulist';
Mumean = mean(Mulist);

%% Theta
Theta1mask = Theta1.*maskx; Theta1list = Theta1(nn>0);
Theta2mask = Theta2.*maskx; Theta2list = Theta2(nn>0);

ThlistAP(1:length(Theta1list),1) = Theta1list;
ThlistAP(1:length(Theta1list),2) = Theta2list;

%% Amplitude weightings
A1p = A1./(A1+A2); A2p = A2./(A1+A2);
A1mask = A1p.*maskx; A1mask(isnan(A1mask)) = 0;
A1list = A1mask(nn>0);
A2mask = A2p.*maskx; A2mask(isnan(A2mask)) = 0;
A2list = A2mask(nn>0);

AmplistAP(1:length(A1list),1) = A1list;
AmplistAP(1:length(A1list),2) = A2list;

%% FA
FAmask = FA.*maskx; FAlist = FA(nn>0);
FAtractlistAP(1:length(FAlist)) = FAlist';

%% Slow/Fast
UsV1list = abs(UsV1(nn>0)); UsV1list = UsV1list(~isnan(UsV1list));
UfV1list = abs(UfV1(nn>0)); UfV1list = UfV1list(~isnan(UfV1list));
UsV2list = abs(UsV2(nn>0)); UsV2list = UsV2list(~isnan(UsV2list));
UfV2list = abs(UfV2(nn>0)); UfV2list = UfV2list(~isnan(UfV2list));

UspercV1 = UsV1list./(UfV1list+UsV1list)*100; UfpercV1 = UfV1list./(UfV1list+UsV1list)*100;
UspercV2 = UsV2list./(UfV2list+UsV2list)*100; UfpercV2 = UfV2list./(UfV2list+UsV2list)*100;

UperctotAP(1:length(UspercV1),3) = UspercV2; UperctotAP(1:length(UspercV1),4) = UfpercV2;
UperctotAP(1:length(UspercV1),1) = UspercV1; UperctotAP(1:length(UspercV1),2) = UfpercV1;

Us = flip(flip(permute(U_s>U_f,[2 1 3]),1),2); Uf = flip(flip(permute(U_f>U_s,[2 1 3]),1),2);
Usmask = maskx.*Us; Ufmask = maskx.*Uf; masktot = Usmask+Ufmask;
Slow = Us(Usmask == 1); Slowperc = length(Slow)/length(maskx(maskx==1))*100;
Fast = Uf(Ufmask == 1); Fastperc = length(Fast)/length(maskx(maskx==1))*100;

TotalAP(1) = Mumean; TotalAP(2) = Slowperc; TotalAP(3) = Fastperc;

%disp(sprintf('Mulist: %i; UsV1list: %i; Theta1list: %i; ii: %i',length(Mulist),length(UsV1list),length(Theta1list)))

save('Quant_MUFS_AP.mat','TotalAP','UperctotAP','MutotlistAP','ThlistAP','FAtractlistAP','AmplistAP')

%% LR
cd ..; cd XBOXXZ
%addpath('PostInv')
load Mu_NLI.mat; load slowfasttwo.mat; load BuckyOut.mat;

%Vmult1 = flip(flip(permute(Vmult,[2 1 3]),1),2);
Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
Theta1 = acos(dot(avec_img,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(avec_img,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));
A1 = flip(flip(permute(A1x,[2 1 3]),1),2); 
A2 = flip(flip(permute(A2x,[2 1 3]),1),2);
UsV1 = flip(flip(permute(U_s_V1,[2 1 3]),1),2); 
UfV1 = flip(flip(permute(U_f_V1,[2 1 3]),1),2);
UsV2 = flip(flip(permute(U_s_V2,[2 1 3]),1),2); 
UfV2 = flip(flip(permute(U_f_V2,[2 1 3]),1),2);

maskx = DTI.*maskXY.*maskXZ;
list = maskx(maskx==1);

%% Mu

nn = Mux.*maskx; nn(isnan(nn))  = 0;
Mulist = nn(nn>0);
figure;im(nn);caxis([0 5000])
MutotlistLR(1:length(Mulist)) = Mulist';
Mumean = mean(Mulist);

%% Theta
Theta1mask = Theta1.*maskx; Theta1list = Theta1(nn>0);
Theta2mask = Theta2.*maskx; Theta2list = Theta2(nn>0);

ThlistLR(1:length(Theta1list),1,ii) = Theta1list;
ThlistLR(1:length(Theta1list),2,ii) = Theta2list;

%% Amplitude weightings
A1p = A1./(A1+A2); A2p = A2./(A1+A2);
A1mask = A1p.*maskx; A1mask(isnan(A1mask)) = 0;
A1list = A1mask(nn>0);
A2mask = A2p.*maskx; A2mask(isnan(A2mask)) = 0;
A2list = A2mask(nn>0);

AmplistLR(1:length(A1list),1) = A1list;
AmplistLR(1:length(A1list),2) = A2list;

%% FA
FAmask = FA.*maskx; FAlist = FA(nn>0);
FAtractlistLR(1:length(FAlist)) = FAlist';

%% Slow/Fast
UsV1list = abs(UsV1(nn>0)); UsV1list = UsV1list(~isnan(UsV1list));
UfV1list = abs(UfV1(nn>0)); UfV1list = UfV1list(~isnan(UfV1list));
UsV2list = abs(UsV2(nn>0)); UsV2list = UsV2list(~isnan(UsV2list));
UfV2list = abs(UfV2(nn>0)); UfV2list = UfV2list(~isnan(UfV2list));

UspercV1 = UsV1list./(UfV1list+UsV1list)*100; UfpercV1 = UfV1list./(UfV1list+UsV1list)*100;
UspercV2 = UsV2list./(UfV2list+UsV2list)*100; UfpercV2 = UfV2list./(UfV2list+UsV2list)*100;

UperctotLR(1:length(UspercV1),3) = UspercV2; UperctotLR(1:length(UspercV1),4) = UfpercV2;
UperctotLR(1:length(UspercV1),1) = UspercV1; UperctotLR(1:length(UspercV1),2) = UfpercV1;

Us = flip(flip(permute(U_s>U_f,[2 1 3]),1),2); Uf = flip(flip(permute(U_f>U_s,[2 1 3]),1),2);
Usmask = maskx.*Us; Ufmask = maskx.*Uf; masktot = Usmask+Ufmask;
Slow = Us(Usmask == 1); Slowperc = length(Slow)/length(maskx(maskx==1))*100;
Fast = Uf(Ufmask == 1); Fastperc = length(Fast)/length(maskx(maskx==1))*100;

TotalLR(1) = Mumean; TotalLR(2) = Slowperc; TotalLR(3) = Fastperc;

save('Quant_MUFS_LR.mat','TotalLR','UperctotLR','MutotlistLR','ThlistLR','FAtractlistLR','AmplistLR')
