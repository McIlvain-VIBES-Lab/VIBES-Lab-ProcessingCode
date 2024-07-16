function quant_perc20_new()
%% Load in AP/LR files from their seperate locations
% Change Locations to AP and LR folders
load Buckyout.mat
load slowfasttwo.mat
dirs = dir('TBI*');
load(dirs.name)
[percent20_AP,Unorms_AP,perc20_AP] = percent(A1x,A2x,U_s_V1,U_f_V1,U_s_V2,U_f_V2,mask);
U_s = U_s_V1+U_s_V2;
totalvox = sum(sum(sum(abs(U_s)>0)));
V1_AP = V1x; V2_AP = V2x;

cd ../LR
load Buckyout.mat
load slowfasttwo.mat
dirs = dir('TBI*');
load(dirs.name)
[percent20_LR,Unorms_LR,perc20_LR] = percent(A1x,A2x,U_s_V1,U_f_V1,U_s_V2,U_f_V2,mask);
V1_LR = V1x; V2_LR = V2x;

cd ../../Tracts/WMmain
CCx = load_nii('tractsCC.nii'); CC = CCx.img>0.25; CC = flip(flip(permute(CC,[2 1 3]),1),2);
CRx = load_nii('tractsCR.nii'); CR = CRx.img>0.25; CR = flip(flip(permute(CR,[2 1 3]),1),2);
CCSx = load_nii('tractsCCS.nii'); CCS = CCSx.img>0.25; CCS = flip(flip(permute(CCS,[2 1 3]),1),2);
SLFx = load_nii('tractsSLF.nii'); SLF = SLFx.img>0.25; SLF = flip(flip(permute(SLF,[2 1 3]),1),2);

cd ../WMminor/niiTracts
CCbx = load_nii('BodyCC.nii'); CCb = CCbx.img>0.25; CCbody = flip(flip(permute(CCb,[2 1 3]),1),2);
cd ..; cd ..; cd ..

%% Theta Differences between Wave propagations
% Theta is the inverse cosine of the dot product of the two vectors 

diff_V1AP_V1LR = rad2deg(abs(acos(dot(V1_AP,V1_LR,4))))>15; % AP1 vs LR1
diff_V1AP_V2AP = rad2deg(abs(acos(dot(V1_AP,V2_AP,4))))>15; % AP1 vs AP2
diff_V1AP_V2LR = rad2deg(abs(acos(dot(V1_AP,V2_LR,4))))>15; % AP1 vs LR2
diff_V1LR_V2AP = rad2deg(abs(acos(dot(V1_LR,V2_AP,4))))>15; % LR1 vs AP2
diff_V1LR_V2LR = rad2deg(abs(acos(dot(V1_LR,V2_LR,4))))>15; % LR1 vs LR2
diff_V2AP_V2LR = rad2deg(abs(acos(dot(V2_AP,V2_LR,4))))>15; % AP2 vs LR2

%% Slow/Fast percentages 

Us1norm_AP = abs(Unorms_AP(:,:,:,1)); Uf1norm_AP = abs(Unorms_AP(:,:,:,2)); Us2norm_AP = abs(Unorms_AP(:,:,:,3)); Uf2norm_AP = abs(Unorms_AP(:,:,:,4));
Us1norm_LR = abs(Unorms_LR(:,:,:,1)); Uf1norm_LR = abs(Unorms_LR(:,:,:,2)); Us2norm_LR = abs(Unorms_LR(:,:,:,3)); Uf2norm_LR = abs(Unorms_LR(:,:,:,4));
Uall_AP = Us1norm_AP+Uf1norm_AP+Us2norm_AP+Uf2norm_AP; Uall_LR = Us1norm_LR+Uf1norm_LR+Us2norm_LR+Uf2norm_LR; Uall_SI = Us1norm_SI+Uf1norm_SI+Us2norm_SI+Uf2norm_SI;

% Slow and fast percentages
Usperc_AP = (Us1norm_AP+Us2norm_AP)./Uall_AP; Ufperc_AP = (Uf1norm_AP+Uf2norm_AP)./Uall_AP;
Usperc_LR = (Us1norm_LR+Us2norm_LR)./Uall_LR; Ufperc_LR = (Uf1norm_LR+Uf2norm_LR)./Uall_LR;

Usperc = cat(4,Usperc_AP,Usperc_LR); Ufperc = cat(4,Ufperc_AP,Ufperc_LR); Uperc = cat(5,Usperc,Ufperc);

%% Slow/Fast waves amounts over 20% of wave magnitude
Us1perc20_AP = perc20_AP(:,:,:,1); Uf1perc20_AP = perc20_AP(:,:,:,2); Us2perc20_AP = perc20_AP(:,:,:,3); Uf2perc20_AP = perc20_AP(:,:,:,4);
Us1perc20_LR = perc20_LR(:,:,:,1); Uf1perc20_LR = perc20_LR(:,:,:,2); Us2perc20_LR = perc20_LR(:,:,:,3); Uf2perc20_LR = perc20_LR(:,:,:,4);

Us1perc_all = Us1perc20_AP+Us1perc20_LR; Uf1perc_all = Uf1perc20_AP+Uf1perc20_LR;
Us2perc_all = Us2perc20_AP+Us2perc20_LR; Uf2perc_all = Uf2perc20_AP+Uf2perc20_LR;

UsAP = Us1perc20_AP+Us2perc20_AP; UfAP = Uf1perc20_AP+Uf2perc20_AP;
UsLR = Us1perc20_LR+Us2perc20_LR; UfLR = Uf1perc20_LR+Uf2perc20_LR;
Uslowall = cat(4,UsAP,UsLR); Ufastall = cat(4,UfAP,UfLR);
UAPall = Us1perc20_AP+Uf1perc20_AP+Us2perc20_AP+Uf2perc20_AP; ULRall = Us1perc20_LR+Uf1perc20_LR+Us2perc20_LR+Uf2perc20_LR;
Usall = Us1perc_all+Us2perc_all; Ufall = Uf1perc_all+Uf2perc_all;

UAPallperc = sum(sum(sum(UAPall==4)))/totalvox; ULRallperc = sum(sum(sum(ULRall==4)))/totalvox;
Uallslowperc = sum(sum(sum(Usall>=2)))/totalvox;Uallfastperc = sum(sum(sum(Ufall>=2)))/totalvox;
Uallallperc = sum(sum(sum((Usall>=2).*(Ufall>=2))))/totalvox;

Uall = [UAPallperc,ULRallperc,Uallallperc];

%% Tweten Criteria for wave magnitudes inclusion (differentiating between waves)
V1s_AP = zeros(80,80,48); V1f_AP = zeros(80,80,48); V2s_AP = zeros(80,80,48); V2f_AP = zeros(80,80,48);
V1s_LR = zeros(80,80,48); V1f_LR = zeros(80,80,48); V2s_LR = zeros(80,80,48); V2f_LR = zeros(80,80,48);

% 20 percent Tweten criteria
for ii = 1:80
    for jj = 1:80
        for kk = 1:48
            if Us1perc20_AP(ii,jj,kk) == 1
                V1s_AP(ii,jj,kk) = 1;
            end
            if Uf1perc20_AP(ii,jj,kk) == 1
                V1f_AP(ii,jj,kk) = 1;
            end
            if Us2perc20_AP(ii,jj,kk) == 1
                V2s_AP(ii,jj,kk) = 1;
            end
            if Uf2perc20_AP(ii,jj,kk) == 1
                V2f_AP(ii,jj,kk) = 1;
            end
            if Us1perc20_LR(ii,jj,kk) == 1 && diff_V1AP_V1LR(ii,jj,kk) == 1 && diff_V1LR_V2AP(ii,jj,kk) == 1
                V1s_LR(ii,jj,kk) = 1;
            end
            if Uf1perc20_LR(ii,jj,kk) == 1 && diff_V1AP_V1LR(ii,jj,kk) == 1 && diff_V1LR_V2AP(ii,jj,kk) == 1
                V1f_LR(ii,jj,kk) = 1;
            end
            if Us2perc20_LR(ii,jj,kk) == 1 && diff_V1AP_V2LR(ii,jj,kk) == 1 && diff_V2AP_V2LR(ii,jj,kk) == 1
                V2s_LR(ii,jj,kk) = 1;
            end
            if Uf2perc20_LR(ii,jj,kk) == 1 && diff_V1AP_V2LR(ii,jj,kk) == 1 && diff_V2AP_V2LR(ii,jj,kk) == 1
                V2f_LR(ii,jj,kk) = 1;
            end
        end
    end
end

% number of voxels fitting criteria/total number of voxels
V1sAP_cp = sum(sum(sum(V1s_AP)))/totalvox; V1fAP_cp = sum(sum(sum(V1f_AP)))/totalvox; V2sAP_cp = sum(sum(sum(V2s_AP)))/totalvox; V2fAP_cp = sum(sum(sum(V2f_AP)))/totalvox;
V1sLR_cp = sum(sum(sum(V1s_LR)))/totalvox; V1fLR_cp = sum(sum(sum(V1f_LR)))/totalvox; V2sLR_cp = sum(sum(sum(V2s_LR)))/totalvox; V2fLR_cp = sum(sum(sum(V2f_LR)))/totalvox;

VtotAP = V1s_AP+V1f_AP+V2s_AP+V2f_AP; VslowAP = V1s_AP+V2s_AP; VfastAP = V1f_AP+V2f_AP;
VtotLR = V1s_LR+V1f_LR+V2s_LR+V2f_LR; VslowLR = V1s_LR+V2s_LR; VfastLR = V1f_LR+V2f_LR;
VslowAP1 = Us1perc20_AP+Us2perc20_AP; VfastAP1 = Uf1perc20_AP+Uf2perc20_AP;
VslowLR1 = Us1perc20_LR+Us2perc20_LR; VfastLR1 = Uf1perc20_LR+Uf2perc20_LR;
Vslowall_c = V1s_AP+V1s_LR+V2s_AP+V2s_LR; Vfastall_c = V1f_AP+V1f_LR+V2f_AP+V2f_LR;
Vall_c = Vslowall_c+Vfastall_c; Vmult = (Vslowall_c>=2).*(Vfastall_c>=2);

VAPall_cp = sum(sum(sum(VtotAP==4)))/totalvox; VLRall_cp = sum(sum(sum(VtotLR==4)))/totalvox;
Vslowperc = sum(sum(sum(Vslowall_c>=2)))/totalvox; Vfastperc = sum(sum(sum(Vfastall_c>=2)))/totalvox;
Vallallperc = sum(sum(sum((Vslowall_c>=2).*(Vfastall_c>=2))))/totalvox;

TW_cps = [VAPall_cp,VLRall_cp,Vallallperc];
Vslow = cat(4,VslowAP,VslowLR,Vslowall_c);
Vfast = cat(4,VfastAP,VfastLR,Vfastall_c);
Vslowtwo_c = V1s_AP+V1s_LR+V2s_AP+V2s_LR; Vfasttwo_c = V1f_AP+V1f_LR+V2f_AP+V2f_LR;

%% Tracts analysis

% All White Matter
Vtotslowall = Vslowtwo_c; Vtotfastall = Vfasttwo_c; totalvoxall = sum(sum(sum(abs(U_s)>0)));
Vslowpercall = sum(sum(sum(Vtotslowall>=2)))/totalvoxall; Vfastpercall = sum(sum(sum(Vtotfastall>=2)))/totalvoxall;
Vallallpercall = sum(sum(sum((Vtotslowall>=2).*(Vtotfastall>=2))))/totalvoxall;
VAPslowall = VslowAP1; VAPfastall = VfastAP1; VLRslowall = VslowLR1; VLRfastall = VfastLR1;
VAPpercall = sum(sum(sum((VAPslowall>=2).*(VAPfastall>=2))))/totalvoxall;
VLRpercall = sum(sum(sum((VLRslowall>=2).*(VLRfastall>=2))))/totalvoxall;

%Body of the Corpus Callosum
VtotslowCCb = Vslowtwo_c.*CCbody; VtotfastCCb = Vfasttwo_c.*CCbody; totalvoxCCb = sum(sum(sum(CCbody.*(abs(U_s)>0))));
VslowpercCCb = sum(sum(sum(VtotslowCCb>=2)))/totalvoxCCb; VfastpercCCb = sum(sum(sum(VtotfastCCb>=2)))/totalvoxCCb;
VallallpercCCb = sum(sum(sum((VtotslowCCb>=2).*(VtotfastCCb>=2))))/totalvoxCCb;
VAPslowCCb = VslowAP1.*CCb; VAPfastCCb = VfastAP1.*CCb; 
VLRslowCCb = VslowLR1.*CCb; VLRfastCCb = VfastLR1.*CCb;
VAPpercCCb = sum(sum(sum((VAPslowCCb>=2).*(VAPfastCCb>=2))))/totalvoxCCb;
VLRpercCCb = sum(sum(sum((VLRslowCCb>=2).*(VLRfastCCb>=2))))/totalvoxCCb;

% Corpus Callosum
VtotslowCC = Vslowtwo_c.*CC; VtotfastCC = Vfasttwo_c.*CC; totalvoxCC = sum(sum(sum(CC.*(abs(U_s)>0))));
VslowpercCC = sum(sum(sum(VtotslowCC>=2)))/totalvoxCC; VfastpercCC = sum(sum(sum(VtotfastCC>=2)))/totalvoxCC;
VallallpercCC = sum(sum(sum((VtotslowCC>=2).*(VtotfastCC>=2))))/totalvoxCC;
VAPslowCC = VslowAP1.*CC; VAPfastCC = VfastAP1.*CC; 
VLRslowCC = VslowLR1.*CC; VLRfastCC = VfastLR1.*CC;
VAPpercCC = sum(sum(sum((VAPslowCC>=2).*(VAPfastCC>=2))))/totalvoxCC;
VLRpercCC = sum(sum(sum((VLRslowCC>=2).*(VLRfastCC>=2))))/totalvoxCC;

% Corona Radiata
VtotslowCR = Vslowtwo_c.*CR; VtotfastCR = Vfasttwo_c.*CR; totalvoxCR = sum(sum(sum(CR.*(abs(U_s)>0))));
VslowpercCR = sum(sum(sum(VtotslowCR>=2)))/totalvoxCR; VfastpercCR = sum(sum(sum(VtotfastCR>=2)))/totalvoxCR;
VallallpercCR = sum(sum(sum((VtotslowCR>=2).*(VtotfastCR>=2))))/totalvoxCR;
VAPslowCR = VslowAP1.*CR; VAPfastCR = VfastAP1.*CR; 
VLRslowCR = VslowLR1.*CR; VLRfastCR = VfastLR1.*CR;
VAPpercCR = sum(sum(sum((VAPslowCR>=2).*(VAPfastCR>=2))))/totalvoxCR;
VLRpercCR = sum(sum(sum((VLRslowCR>=2).*(VLRfastCR>=2))))/totalvoxCR;

% Corticolspinal Tract
VtotslowCCS = Vslowtwo_c.*CCS; VtotfastCCS = Vfasttwo_c.*CCS; totalvoxCCS = sum(sum(sum(CCS.*(abs(U_s)>0))));
VslowpercCCS = sum(sum(sum(VtotslowCCS>=2)))/totalvoxCCS; VfastpercCCS = sum(sum(sum(VtotfastCCS>=2)))/totalvoxCCS;
VallallpercCCS = sum(sum(sum((VtotslowCCS>=2).*(VtotfastCCS>=2))))/totalvoxCCS;
VAPslowCCS = VslowAP1.*CCS; VAPfastCCS = VfastAP1.*CCS; 
VLRslowCCS = VslowLR1.*CCS; VLRfastCCS = VfastLR1.*CCS;
VAPpercCCS = sum(sum(sum((VAPslowCCS>=2).*(VAPfastCCS>=2))))/totalvoxCCS;
VLRpercCCS = sum(sum(sum((VLRslowCCS>=2).*(VLRfastCCS>=2))))/totalvoxCCS;

% Superior Longitudinal Fasciculus
VtotslowSLF = Vslowtwo_c.*SLF; VtotfastSLF = Vfasttwo_c.*SLF; totalvoxSLF = sum(sum(sum(SLF.*(abs(U_s)>0))));
VslowpercSLF = sum(sum(sum(VtotslowSLF>=2)))/totalvoxSLF; VfastpercSLF = sum(sum(sum(VtotfastSLF>=2)))/totalvoxSLF;
VallallpercSLF = sum(sum(sum((VtotslowSLF>=2).*(VtotfastSLF>=2))))/totalvoxSLF;
VAPslowSLF = VslowAP1.*SLF; VAPfastSLF = VfastAP1.*SLF; 
VLRslowSLF = VslowLR1.*SLF; VLRfastSLF = VfastLR1.*SLF;
VAPpercSLF = sum(sum(sum((VAPslowSLF>=2).*(VAPfastSLF>=2))))/totalvoxSLF;
VLRpercSLF = sum(sum(sum((VLRslowSLF>=2).*(VLRfastSLF>=2))))/totalvoxSLF;

TW_tracts = [VallallpercCC,VallallpercCR,VallallpercCCS,VallallpercSLF];
all_all = [VAPpercall,VLRpercall,Vallallpercall];
CCb_all = [VAPpercCCb,VLRpercCCb,VallallpercCCb];
CC_all = [VAPpercCC,VLRpercCC,VallallpercCC];
CR_all = [VAPpercCR,VLRpercCR,VallallpercCR];
CCS_all = [VAPpercCCS,VLRpercCCS,VallallpercCCS];
SLF_all = [VAPpercSLF,VLRpercSLF,VallallpercSLF];

save('SlowFastSummation.mat','all_all','Uall','TW_cps','TW_tracts','CCb_all','CC_all','CR_all','CCS_all','SLF_all','Vslow','Vfast','Uslowall','Ufastall','Vmult')
save('Criteria.mat','V1s_AP','V1f_AP','V2s_AP','V2f_AP','V1s_LR','V1f_LR','V2s_LR','V2f_LR')