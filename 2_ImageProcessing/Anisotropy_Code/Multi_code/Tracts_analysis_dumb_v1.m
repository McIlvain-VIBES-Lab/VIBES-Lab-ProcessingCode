load('dtiall.mat');
%% Load in the Tracts
cd Tracts
CC_ = load_nii('tractsCC.nii'); CC = repmat(CC_.img>.6,[1 1 1 3]); CCT = CC_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CR_ = load_nii('tractsCR.nii'); CR = repmat(CR_.img>.6,[1 1 1 3]); CRT = CR_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CST_ = load_nii('tractsCCS.nii'); CST = repmat(CST_.img>.6, [1 1 1 3]); CSTT = CST_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
SLF_ = load_nii('tractsSLF.nii'); SLF = repmat(SLF_.img>.6, [1 1 1 3]); SLFT = SLF_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
cd ..

%% Load in AP F/S (percentage and Theta dumbed down version)
% Means and Standard Devs
cd MRE; cd AP;
% Load Dumbed down version
load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..
MSvox_AP = MSvox; MFvox_AP = MFvox; 
CC_Svox_AP = MSvox.*CCT; CR_Svox_AP = MSvox.*CRT; CST_Svox_AP = MSvox.*CSTT; SLF_Svox_AP = MSvox.*SLFT;
CC_Fvox_AP = MFvox.*CCT; CR_Fvox_AP = MFvox.*CRT; CST_Fvox_AP = MFvox.*CSTT; SLF_Fvox_AP = MFvox.*SLFT;
CC_Sperc_AP = sum(CC_Svox_AP(CC_Svox_AP==1))/sum(CCT(CCT==1)).*100; CR_Sperc_AP = sum(CR_Svox_AP(CR_Svox_AP==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_AP = sum(CST_Svox_AP(CST_Svox_AP==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_AP = sum(SLF_Svox_AP(SLF_Svox_AP==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_AP = sum(CC_Fvox_AP(CC_Fvox_AP==1))/sum(CCT(CCT==1)).*100; CR_Fperc_AP = sum(CR_Fvox_AP(CR_Fvox_AP==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_AP = sum(CST_Fvox_AP(CST_Fvox_AP==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_AP = sum(SLF_Fvox_AP(SLF_Fvox_AP==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.2
THETA = radtodeg(THETA);
CC_Theta_AP = THETA.*CCT; CR_Theta_AP = THETA.*CRT; CST_Theta_AP = THETA.*CSTT; SLF_Theta_AP = THETA.*SLFT;
CC_Theta_AP_mean = mean(CC_Theta_AP(CC_Theta_AP>0)); CC_Theta_AP_stdev = std(CC_Theta_AP(CC_Theta_AP>0),'omitnan');
CR_Theta_AP_mean = mean(CR_Theta_AP(CR_Theta_AP>0)); CR_Theta_AP_stdev = std(CR_Theta_AP(CR_Theta_AP>0),'omitnan');
CST_Theta_AP_mean = mean(CST_Theta_AP(CST_Theta_AP>0)); CST_Theta_AP_stdev = std(CST_Theta_AP(CST_Theta_AP>0),'omitnan');
SLF_Theta_AP_mean = mean(SLF_Theta_AP(SLF_Theta_AP>0)); SLF_Theta_AP_stdev = std(SLF_Theta_AP(SLF_Theta_AP>0),'omitnan');

%% Load in LR F/S (percentage and Theta dumbed down version)
% Means and Standard Devs
cd MRE; cd LR;
% Load Dumbed down version
load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..
MSvox_LR = MSvox; MFvox_LR = MFvox;
CC_Svox_LR = MSvox.*CCT; CR_Svox_LR = MSvox.*CRT; CST_Svox_LR = MSvox.*CSTT; SLF_Svox_LR = MSvox.*SLFT;
CC_Fvox_LR = MFvox.*CCT; CR_Fvox_LR = MFvox.*CRT; CST_Fvox_LR = MFvox.*CSTT; SLF_Fvox_LR = MFvox.*SLFT;
CC_Sperc_LR = sum(CC_Svox_LR(CC_Svox_LR==1))/sum(CCT(CCT==1)).*100; CR_Sperc_LR = sum(CR_Svox_LR(CR_Svox_LR==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_LR = sum(CST_Svox_LR(CST_Svox_LR==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_LR = sum(SLF_Svox_LR(SLF_Svox_LR==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_LR = sum(CC_Fvox_LR(CC_Fvox_LR==1))/sum(CCT(CCT==1)).*100; CR_Fperc_LR = sum(CR_Fvox_LR(CR_Fvox_LR==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_LR = sum(CST_Fvox_LR(CST_Fvox_LR==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_LR = sum(SLF_Fvox_LR(SLF_Fvox_LR==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.2
THETA = radtodeg(THETA);
CC_Theta_LR = THETA.*CCT; CR_Theta_LR = THETA.*CRT; CST_Theta_LR = THETA.*CSTT; SLF_Theta_LR = THETA.*SLFT;
CC_Theta_LR_mean = mean(CC_Theta_LR(CC_Theta_LR>0)); CC_Theta_LR_stdev = std(CC_Theta_LR(CC_Theta_LR>0),'omitnan');
CR_Theta_LR_mean = mean(CR_Theta_LR(CR_Theta_LR>0)); CR_Theta_LR_stdev = std(CR_Theta_LR(CR_Theta_LR>0),'omitnan');
CST_Theta_LR_mean = mean(CST_Theta_LR(CST_Theta_LR>0)); CST_Theta_LR_stdev = std(CST_Theta_LR(CST_Theta_LR>0),'omitnan');
SLF_Theta_LR_mean = mean(SLF_Theta_LR(SLF_Theta_LR>0)); SLF_Theta_LR_stdev = std(SLF_Theta_LR(SLF_Theta_LR>0),'omitnan');

%% Load in SI F/S (percentage and Theta dumbed down version)
% Means and Standard Devs
cd MRE; cd SI;
% Load Dumbed down version
load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..

MSvox_SI = MSvox; MFvox_SI = MFvox;

CC_Svox_SI = MSvox.*CCT; CR_Svox_SI = MSvox.*CRT; CST_Svox_SI = MSvox.*CSTT; SLF_Svox_SI = MSvox.*SLFT;
CC_Fvox_SI = MFvox.*CCT; CR_Fvox_SI = MFvox.*CRT; CST_Fvox_SI = MFvox.*CSTT; SLF_Fvox_SI = MFvox.*SLFT;
CC_Sperc_SI = sum(CC_Svox_SI(CC_Svox_SI==1))/sum(CCT(CCT==1)).*100; CR_Sperc_SI = sum(CR_Svox_SI(CR_Svox_SI==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_SI = sum(CST_Svox_SI(CST_Svox_SI==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_SI = sum(SLF_Svox_SI(SLF_Svox_SI==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_SI = sum(CC_Fvox_SI(CC_Fvox_SI==1))/sum(CCT(CCT==1)).*100; CR_Fperc_SI = sum(CR_Fvox_SI(CR_Fvox_SI==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_SI = sum(CST_Fvox_SI(CST_Fvox_SI==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_SI = sum(SLF_Fvox_SI(SLF_Fvox_SI==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.2
THETA = radtodeg(THETA);
CC_Theta_SI = THETA.*CCT; CR_Theta_SI = THETA.*CRT; CST_Theta_SI = THETA.*CSTT; SLF_Theta_SI = THETA.*SLFT;
CC_Theta_SI_mean = mean(CC_Theta_SI(CC_Theta_SI>0)); CC_Theta_SI_stdev = std(CC_Theta_SI(CC_Theta_SI>0),'omitnan');
CR_Theta_SI_mean = mean(CR_Theta_SI(CR_Theta_SI>0)); CR_Theta_SI_stdev = std(CR_Theta_SI(CR_Theta_SI>0),'omitnan');
CST_Theta_SI_mean = mean(CST_Theta_SI(CST_Theta_SI>0)); CST_Theta_SI_stdev = std(CST_Theta_SI(CST_Theta_SI>0),'omitnan');
SLF_Theta_SI_mean = mean(SLF_Theta_SI(SLF_Theta_SI>0)); SLF_Theta_SI_stdev = std(SLF_Theta_SI(SLF_Theta_SI>0),'omitnan');

%% Output 

AP_THETA = [CC_Theta_AP_mean CR_Theta_AP_mean CST_Theta_AP_mean SLF_Theta_AP_mean ; CC_Theta_AP_stdev CR_Theta_AP_stdev CST_Theta_AP_stdev SLF_Theta_AP_stdev]
AP_SPERC = [CC_Sperc_AP CR_Sperc_AP CST_Sperc_AP SLF_Sperc_AP]
AP_FPERC = [CC_Fperc_AP CR_Fperc_AP CST_Fperc_AP SLF_Fperc_AP]

LR_THETA = [CC_Theta_LR_mean CR_Theta_LR_mean CST_Theta_LR_mean SLF_Theta_LR_mean ; CC_Theta_LR_stdev CR_Theta_LR_stdev CST_Theta_LR_stdev SLF_Theta_LR_stdev]
LR_SPERC = [CC_Sperc_LR CR_Sperc_LR CST_Sperc_LR SLF_Sperc_LR]
LR_FPERC = [CC_Fperc_LR CR_Fperc_LR CST_Fperc_LR SLF_Fperc_LR]

SI_THETA = [CC_Theta_SI_mean CR_Theta_SI_mean CST_Theta_SI_mean SLF_Theta_SI_mean ; CC_Theta_SI_stdev CR_Theta_SI_stdev CST_Theta_SI_stdev SLF_Theta_SI_stdev]
SI_SPERC = [CC_Sperc_SI CR_Sperc_SI CST_Sperc_SI SLF_Sperc_SI]
SI_FPERC = [CC_Fperc_SI CR_Fperc_SI CST_Fperc_SI SLF_Fperc_SI]