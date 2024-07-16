%% TRACTS
% 1 = Anterior Corona Radiata Left
% 2 = Anterior Corona Radiata Right
% 3 = Body Corpus Callosum
% 4 = Corticolspinal Tract Left
% 5 = Corticolspinal Tract Right
% 6 = Genu Corpus Callosum
% 7 = Posterior Corona Radiata Left
% 8 = Posterior Corona Radiata Right
% 9 = Superior Corona Radiata Left
% 10 = Superior Corona Radiata Right
% 11 = Superior Longitudinal Fasciculus Left
% 12 = Superior Longitudinal Fasciculus Right
% 13 = Splenium Corpus Callosum

%% CC Tracts ç

cd Tracts; cd WMminor; load Tracts.mat
cd niiTracts; niis = dir('*.nii');

CC = double(Tracts(:,:,:,[3 6 13])); 
CC(:,:,:,4) = double(sum(Tracts(:,:,:,[3 6 13]),4)>0);

cd ..; cd ..; cd .. 

load dtiall.mat
load Criteria.mat
DTI = abs(dtiall(:,:,:,1))>0;
% cd Buckley_Tbi_Brain
cd MRE; cd AP
load Multiexcite_03_AP_170222.mat
maskAP = flip(flip(permute(mask,[2 1 3]),1),2);
cd ..; cd LR
load Multiexcite_03_LR_170222.mat
maskLR = flip(flip(permute(mask,[2 1 3]),1),2);
cd ..; cd .. 

%%  AP
% cd Skeleton
load FA.mat
% cd ..
load SlowFastSummation.mat; load Criteria.mat
cd MRE; cd AP 
%addpath('PostInv')
load Mu.mat; load slowfasttwo.mat; load BuckyOut.mat; 
% cd LDI_compiled
load Mu_LDI_diff732.mat

MutotlistAP = zeros(1,4,4); TotalAP = zeros(4,3); AmplistAP = zeros(1,4,2);
UperctotAP = zeros(1,4,4); ThlistAP = zeros(1,2,4); FAtractlistAP = zeros(1,4);
critAP = zeros(1,4,4); GplistAP = zeros(1,4,4);

Vmult1 = flip(flip(permute(Vmult,[2 1 3]),1),2);
Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
Mus1 = Mu_s1; Muf1 = Mu_f1; Mus2 = Mu_s2; Muf2 = Mu_f2;
Theta1 = acos(dot(dtiall,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(dtiall,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));

for ii = 1:4
    vox = zeros(size(Mus1));
    TR = CC(:,:,:,ii); 
    maskx = Vmult1.*TR.*DTI.*maskAP.*maskLR;
    list = maskx(maskx==1);
    
    %% Mu
    nn = Mus1.*maskx; nn(isnan(nn)) = 0;
    Mus1x = Mus1.*maskx; Mus2x = Mus2.*maskx;
    Muf1x = Muf1.*maskx; Muf2x = Muf2.*maskx;
    Mulist = nn(nn>0);
    Mus1list = Mus1x(Mus1x>0); Mus2list = Mus2x(Mus2x>0);
    Muf1list = Muf1x(Muf1x>0); Muf2list = Muf2x(Muf2x>0);
    % figure;im(nn);caxis([0 5000])
    MutotlistAP(1:length(Mus1list),ii,1) = Mus1list;
    MutotlistAP(1:length(Mus2list),ii,2) = Mus2list;
    MutotlistAP(1:length(Muf1list),ii,3) = Muf1list;
    MutotlistAP(1:length(Muf2list),ii,4) = Muf2list;
    
    %% G'
    Gps1 = Gps1x.*maskx; Gps2 = Gps2x.*maskx;
    Gpf1 = Gpf1x.*maskx; Gpf2 = Gpf2x.*maskx;
    Gps1list = Gps1(Gps1>0); Gps2list = Gps2(Gps2>0);
    Gpf1list = Gpf1(Gpf1>0); Gpf2list = Gpf2(Gpf2>0);
    % figure;im(nn);caxis([0 5000])
    GplistAP(1:length(Gps1list),ii,1) = Gps1list;
    GplistAP(1:length(Gps2list),ii,2) = Gps2list;
    GplistAP(1:length(Gpf1list),ii,3) = Gpf1list;
    GplistAP(1:length(Gpf2list),ii,4) = Gpf2list;
    
    %% Theta
    Theta1mask = Theta1.*maskx; Theta1list = Theta1(nn>0);
    Theta2mask = Theta2.*maskx; Theta2list = Theta2(nn>0);
    
    ThlistAP(1:length(Theta1list),1,ii) = Theta1list;
    ThlistAP(1:length(Theta1list),2,ii) = Theta2list;
      
    %% FA
    FAmask = FA.*maskx; FAlist = FA(nn>0);
    FAtractlistAP(1:length(FAlist),ii) = FAlist;
    
    % TotalAP(ii,1) = Mumean; TotalAP(ii,2) = Slowperc; TotalAP(ii,3) = Fastperc;
    V1s = ffp(V1s_AP).*maskx; V1slist = V1s(nn>0);
    V1f = ffp(V1f_AP).*maskx; V1flist = V1f(nn>0);
    V2s = ffp(V2s_AP).*maskx; V2slist = V2s(nn>0);
    V2f = ffp(V2f_AP).*maskx; V2flist = V2f(nn>0);
    critAP(1:length(V1slist),1,ii) = V1slist;
    critAP(1:length(V1flist),2,ii) = V1flist;
    critAP(1:length(V2slist),3,ii) = V2slist;
    critAP(1:length(V2flist),4,ii) = V2flist;
    
    disp(sprintf('Mulist: %i; Theta1list: %i; ii: %i',length(Mulist),length(Theta1list),ii))
    
end

save('Quant_MUFS_APcrit.mat','TotalAP','MutotlistAP','ThlistAP','FAtractlistAP','critAP','GplistAP')
% cd ..
%% LR
cd ..; cd LR
%addpath('PostInv')
load slowfasttwo.mat; load BuckyOut.mat; 
load Mu.mat
% cd LDI_compiled
load Mu_LDI_diff732.mat
MutotlistLR = zeros(1,4,4); TotalLR = zeros(4,3); ThlistLR = zeros(1,2,4); FAtractlistLR = zeros(1,4);
critLR = zeros(1,4,4); GplistLR = zeros(1,4,4);

Vmult1 = flip(flip(permute(Vmult,[2 1 3]),1),2);
Mux = flip(flip(permute(Mu,[2 1 3]),1),2);
Mus1 = Mu_s1; Muf1 = Mu_f1; Mus2 = Mu_s2; Muf2 = Mu_f2;
Theta1 = acos(dot(dtiall,flip(flip(permute(V1x,[2 1 3 4]),1),2),4));
Theta2 = acos(dot(dtiall,flip(flip(permute(V2x,[2 1 3 4]),1),2),4));

for ii = 1:4
    vox = zeros(size(Mus1));
    TR = CC(:,:,:,ii); 
    maskx = Vmult1.*TR.*DTI.*maskAP.*maskLR;
    list = maskx(maskx==1);
    
    %% Mu
    nn = Mus1.*maskx; nn(isnan(nn)) = 0;
    Mus1x = Mus1.*maskx; Mus2x = Mus2.*maskx;
    Muf1x = Muf1.*maskx; Muf2x = Muf2.*maskx;
    Mulist = nn(nn>0);
    Mus1list = Mus1x(Mus1x>0); Mus2list = Mus2x(Mus2x>0);
    Muf1list = Muf1x(Muf1x>0); Muf2list = Muf2x(Muf2x>0);
    % figure;im(nn);caxis([0 5000])
    MutotlistLR(1:length(Mus1list),ii,1) = Mus1list;
    MutotlistLR(1:length(Mus2list),ii,2) = Mus2list;
    MutotlistLR(1:length(Muf1list),ii,3) = Muf1list;
    MutotlistLR(1:length(Muf2list),ii,4) = Muf2list;
    
    
    %% G'
    Gps1 = Gps1x.*maskx; Gps2 = Gps2x.*maskx;
    Gpf1 = Gpf1x.*maskx; Gpf2 = Gpf2x.*maskx;
    Gps1list = Gps1(Gps1>0); Gps2list = Gps2(Gps2>0);
    Gpf1list = Gpf1(Gpf1>0); Gpf2list = Gpf2(Gpf2>0);
    % figure;im(nn);caxis([0 5000])
    GplistLR(1:length(Gps1list),ii,1) = Gps1list;
    GplistLR(1:length(Gps2list),ii,2) = Gps2list;
    GplistLR(1:length(Gpf1list),ii,3) = Gpf1list;
    GplistLR(1:length(Gpf2list),ii,4) = Gpf2list;
    
    
    %% Theta
    Theta1mask = Theta1.*maskx; Theta1list = Theta1(nn>0);
    Theta2mask = Theta2.*maskx; Theta2list = Theta2(nn>0);
    
    ThlistLR(1:length(Theta1list),1,ii) = Theta1list;
    ThlistLR(1:length(Theta1list),2,ii) = Theta2list;
    
    %% FA
    FAmask = FA.*maskx; FAlist = FA(nn>0);
    FAtractlistLR(1:length(FAlist),ii) = FAlist;
       
    % TotalLR(ii,1) = Mumean; TotalLR(ii,2) = Slowperc; TotalLR(ii,3) = Fastperc;
    
    V1s = ffp(V1s_LR).*maskx; V1slist = V1s(nn>0);
    V1f = ffp(V1f_LR).*maskx; V1flist = V1f(nn>0);
    V2s = ffp(V2s_LR).*maskx; V2slist = V2s(nn>0);
    V2f = ffp(V2f_LR).*maskx; V2flist = V2f(nn>0);
    critLR(1:length(V1slist),1,ii) = V1slist;
    critLR(1:length(V1flist),2,ii) = V1flist;
    critLR(1:length(V2slist),3,ii) = V2slist;
    critLR(1:length(V2flist),4,ii) = V2flist;
    
    disp(sprintf('Mulist: %i; Theta1list: %i; ii: %i',length(Mulist),length(Theta1list),ii))
    
end

save('Quant_MUFS_LRcrit.mat','TotalLR','MutotlistLR','ThlistLR','FAtractlistLR','critLR','GplistLR')

cd ..; cd ..; cd ..