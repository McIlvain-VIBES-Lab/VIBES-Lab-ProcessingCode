load('dtiall.mat');
%% Load in the Tracts
cd Tracts
CC_ = load_nii('tractsCC.nii'); CC = repmat(CC_.img>.6,[1 1 1 3]); CCT = CC_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CR_ = load_nii('tractsCR.nii'); CR = repmat(CR_.img>.6,[1 1 1 3]); CRT = CR_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CST_ = load_nii('tractsCCS.nii'); CST = repmat(CST_.img>.6, [1 1 1 3]); CSTT = CST_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
SLF_ = load_nii('tractsSLF.nii'); SLF = repmat(SLF_.img>.6, [1 1 1 3]); SLFT = SLF_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
cd ..

%% Load in AP SNR, VSS, and F/S (percentage and Theta)
% Means and Standard Devs
cd MRE; cd AP; cd PreInv; cd mre_process; cd dcfiles
% Load SNR
load OSS_SNR.mat; cd ..; cd ..; cd ..
AP_SNR = OSS_SNR;
% Load VSS
cd PostInv; load Mu.mat; cd ..
% Load Iterative Data (Tracts percentage) and Predefined data (Theta)
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd Predefined_Outputs; load Outputs.mat; cd .. % Theta
% Load Dumbed down version
%load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..

%% Calculate AP VSS and F/S identity based upon Tracts
% VSS
CC_VSS_AP = Mu.*CCT; CR_VSS_AP = Mu.*CRT; CST_VSS_AP = Mu.*CSTT; SLF_VSS_AP = Mu.*SLFT;

CC_VSS_AP_mean = mean(CC_VSS_AP(CC_VSS_AP>0)); CC_VSS_AP_stdev = std(CC_VSS_AP(CC_VSS_AP>0),'omitnan');
CR_VSS_AP_mean = mean(CR_VSS_AP(CR_VSS_AP>0)); CR_VSS_AP_stdev = std(CR_VSS_AP(CR_VSS_AP>0),'omitnan');
CST_VSS_AP_mean = mean(CST_VSS_AP(CST_VSS_AP>0)); CST_VSS_AP_stdev = std(CST_VSS_AP(CST_VSS_AP>0),'omitnan');
SLF_VSS_AP_mean = mean(SLF_VSS_AP(SLF_VSS_AP>0)); SLF_VSS_AP_stdev = std(SLF_VSS_AP(SLF_VSS_AP>0),'omitnan');

%F/S percentage p.1
CC_Svox_AP = S_vox_U.*CCT; CR_Svox_AP = S_vox_U.*CRT; CST_Svox_AP = S_vox_U.*CSTT; SLF_Svox_AP = S_vox_U.*SLFT;
CC_Fvox_AP = F_vox_U.*CCT; CR_Fvox_AP = F_vox_U.*CRT; CST_Fvox_AP = F_vox_U.*CSTT; SLF_Fvox_AP = F_vox_U.*SLFT;
CC_Sperc_AP = sum(CC_Svox_AP(CC_Svox_AP==1))/sum(CCT(CCT==1)).*100; CR_Sperc_AP = sum(CR_Svox_AP(CR_Svox_AP==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_AP = sum(CST_Svox_AP(CST_Svox_AP==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_AP = sum(SLF_Svox_AP(SLF_Svox_AP==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_AP = sum(CC_Fvox_AP(CC_Fvox_AP==1))/sum(CCT(CCT==1)).*100; CR_Fperc_AP = sum(CR_Fvox_AP(CR_Fvox_AP==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_AP = sum(CST_Fvox_AP(CST_Fvox_AP==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_AP = sum(SLF_Fvox_AP(SLF_Fvox_AP==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.1
CC_Theta_AP = Theta.*CCT; CR_Theta_AP = Theta.*CRT; CST_Theta_AP = Theta.*CSTT; SLF_Theta_AP = Theta.*SLFT;
CC_Theta_AP_mean = mean(CC_Theta_AP(CC_Theta_AP>0)); CC_Theta_AP_stdev = std(CC_Theta_AP(CC_Theta_AP>0),'omitnan');
CR_Theta_AP_mean = mean(CR_Theta_AP(CR_Theta_AP>0)); CR_Theta_AP_stdev = std(CR_Theta_AP(CR_Theta_AP>0),'omitnan');
CST_Theta_AP_mean = mean(CST_Theta_AP(CST_Theta_AP>0)); CST_Theta_AP_stdev = std(CST_Theta_AP(CST_Theta_AP>0),'omitnan');
SLF_Theta_AP_mean = mean(SLF_Theta_AP(SLF_Theta_AP>0)); SLF_Theta_AP_stdev = std(SLF_Theta_AP(SLF_Theta_AP>0),'omitnan');

%% Load in LR SNR, VSS, and tracts (percentage and Theta)
% Means and Standard Devs
cd MRE; cd LR; cd PreInv; cd mre_process; cd dcfiles
% Load SNR
load OSS_SNR.mat; cd ..; cd ..; cd ..
LR_SNR = OSS_SNR;
% Load VSS
cd PostInv; load Mu.mat; cd ..
% Load Iterative Data (Tracts percentage) and Predefined data (Theta)
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd Predefined_Outputs; load Outputs.mat; cd .. % Theta
% Load Dumbed down version
%load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..

%% Calculate LR VSS and F/S identity based upon Tracts

CC_VSS_LR = Mu.*CCT; CR_VSS_LR = Mu.*CRT; CST_VSS_LR = Mu.*CSTT; SLF_VSS_LR = Mu.*SLFT;

CC_VSS_LR_mean = mean(CC_VSS_LR(CC_VSS_LR>0)); CC_VSS_LR_stdev = std(CC_VSS_LR(CC_VSS_LR>0),'omitnan');
CR_VSS_LR_mean = mean(CR_VSS_LR(CR_VSS_LR>0)); CR_VSS_LR_stdev = std(CR_VSS_LR(CR_VSS_LR>0),'omitnan');
CST_VSS_LR_mean = mean(CST_VSS_LR(CST_VSS_LR>0)); CST_VSS_LR_stdev = std(CST_VSS_LR(CST_VSS_LR>0),'omitnan');
SLF_VSS_LR_mean = mean(SLF_VSS_LR(SLF_VSS_LR>0)); SLF_VSS_LR_stdev = std(SLF_VSS_LR(SLF_VSS_LR>0),'omitnan');

%F/S percentage p.1
CC_Svox_LR = S_vox_U.*CCT; CR_Svox_LR = S_vox_U.*CRT; CST_Svox_LR = S_vox_U.*CSTT; SLF_Svox_LR = S_vox_U.*SLFT;
CC_Fvox_LR = F_vox_U.*CCT; CR_Fvox_LR = F_vox_U.*CRT; CST_Fvox_LR = F_vox_U.*CSTT; SLF_Fvox_LR = F_vox_U.*SLFT;
CC_Sperc_LR = sum(CC_Svox_LR(CC_Svox_LR==1))/sum(CCT(CCT==1)).*100; CR_Sperc_LR = sum(CR_Svox_LR(CR_Svox_LR==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_LR = sum(CST_Svox_LR(CST_Svox_LR==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_LR = sum(SLF_Svox_LR(SLF_Svox_LR==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_LR = sum(CC_Fvox_LR(CC_Fvox_LR==1))/sum(CCT(CCT==1)).*100; CR_Fperc_LR = sum(CR_Fvox_LR(CR_Fvox_LR==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_LR = sum(CST_Fvox_LR(CST_Fvox_LR==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_LR = sum(SLF_Fvox_LR(SLF_Fvox_LR==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.1
CC_Theta_LR = Theta.*CCT; CR_Theta_LR = Theta.*CRT; CST_Theta_LR = Theta.*CSTT; SLF_Theta_LR = Theta.*SLFT;
CC_Theta_LR_mean = mean(CC_Theta_LR(CC_Theta_LR>0)); CC_Theta_LR_stdev = std(CC_Theta_LR(CC_Theta_LR>0),'omitnan');
CR_Theta_LR_mean = mean(CR_Theta_LR(CR_Theta_LR>0)); CR_Theta_LR_stdev = std(CR_Theta_LR(CR_Theta_LR>0),'omitnan');
CST_Theta_LR_mean = mean(CST_Theta_LR(CST_Theta_LR>0)); CST_Theta_LR_stdev = std(CST_Theta_LR(CST_Theta_LR>0),'omitnan');
SLF_Theta_LR_mean = mean(SLF_Theta_LR(SLF_Theta_LR>0)); SLF_Theta_LR_stdev = std(SLF_Theta_LR(SLF_Theta_LR>0),'omitnan');


%% Load in SI SNR, VSS, and tracts (percentage and Theta)
% Means and Standard Devs
cd MRE; cd SI; cd PreInv; cd mre_process; cd dcfiles
% Load SNR
load OSS_SNR.mat; cd ..; cd ..; cd ..
SI_SNR = OSS_SNR;
% Load VSS
cd PostInv; load Mu.mat; cd ..
% Load Iterative Data (Tracts percentage) and Predefined data (Theta)
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd Predefined_Outputs; load Outputs.mat; cd .. % Theta
% Load Dumbed down version
%load WashU_dumb.mat %MSLOW, MFAST, THETA, MSvox, MFvox
cd ..; cd ..

%% Calculate SI VSS and F/S identity based upon Tracts

CC_VSS_SI = Mu.*CCT; CR_VSS_SI = Mu.*CRT; CST_VSS_SI = Mu.*CSTT; SLF_VSS_SI = Mu.*SLFT;

CC_VSS_SI_mean = mean(CC_VSS_SI(CC_VSS_SI>0)); CC_VSS_SI_stdev = std(CC_VSS_SI(CC_VSS_SI>0),'omitnan');
CR_VSS_SI_mean = mean(CR_VSS_SI(CR_VSS_SI>0)); CR_VSS_SI_stdev = std(CR_VSS_SI(CR_VSS_SI>0),'omitnan');
CST_VSS_SI_mean = mean(CST_VSS_SI(CST_VSS_SI>0)); CST_VSS_SI_stdev = std(CST_VSS_SI(CST_VSS_SI>0),'omitnan');
SLF_VSS_SI_mean = mean(SLF_VSS_SI(SLF_VSS_SI>0)); SLF_VSS_SI_stdev = std(SLF_VSS_SI(SLF_VSS_SI>0),'omitnan');

%F/S percentage p.1
CC_Svox_SI = S_vox_U.*CCT; CR_Svox_SI = S_vox_U.*CRT; CST_Svox_SI = S_vox_U.*CSTT; SLF_Svox_SI = S_vox_U.*SLFT;
CC_Fvox_SI = F_vox_U.*CCT; CR_Fvox_SI = F_vox_U.*CRT; CST_Fvox_SI = F_vox_U.*CSTT; SLF_Fvox_SI = F_vox_U.*SLFT;
CC_Sperc_SI = sum(CC_Svox_SI(CC_Svox_SI==1))/sum(CCT(CCT==1)).*100; CR_Sperc_SI = sum(CR_Svox_SI(CR_Svox_SI==1))/sum(CRT(CRT==1)).*100;
CST_Sperc_SI = sum(CST_Svox_SI(CST_Svox_SI==1))/sum(CSTT(CSTT==1)).*100; SLF_Sperc_SI = sum(SLF_Svox_SI(SLF_Svox_SI==1))/sum(SLFT(SLFT==1)).*100;
CC_Fperc_SI = sum(CC_Fvox_SI(CC_Fvox_SI==1))/sum(CCT(CCT==1)).*100; CR_Fperc_SI = sum(CR_Fvox_SI(CR_Fvox_SI==1))/sum(CRT(CRT==1)).*100;
CST_Fperc_SI = sum(CST_Fvox_SI(CST_Fvox_SI==1))/sum(CSTT(CSTT==1)).*100; SLF_Fperc_SI = sum(SLF_Fvox_SI(SLF_Fvox_SI==1))/sum(SLFT(SLFT==1)).*100;

%F/S Theta p.1
CC_Theta_SI = Theta.*CCT; CR_Theta_SI = Theta.*CRT; CST_Theta_SI = Theta.*CSTT; SLF_Theta_SI = Theta.*SLFT;
CC_Theta_SI_mean = mean(CC_Theta_SI(CC_Theta_SI>0)); CC_Theta_SI_stdev = std(CC_Theta_SI(CC_Theta_SI>0),'omitnan');
CR_Theta_SI_mean = mean(CR_Theta_SI(CR_Theta_SI>0)); CR_Theta_SI_stdev = std(CR_Theta_SI(CR_Theta_SI>0),'omitnan');
CST_Theta_SI_mean = mean(CST_Theta_SI(CST_Theta_SI>0)); CST_Theta_SI_stdev = std(CST_Theta_SI(CST_Theta_SI>0),'omitnan');
SLF_Theta_SI_mean = mean(SLF_Theta_SI(SLF_Theta_SI>0)); SLF_Theta_SI_stdev = std(SLF_Theta_SI(SLF_Theta_SI>0),'omitnan');

%% Output

SNR = [AP_SNR LR_SNR]

AP_VSS = [CC_VSS_AP_mean CR_VSS_AP_mean CST_VSS_AP_mean SLF_VSS_AP_mean ; CC_VSS_AP_stdev CR_VSS_AP_stdev CST_VSS_AP_stdev SLF_VSS_AP_stdev]
AP_THETA = [CC_Theta_AP_mean CR_Theta_AP_mean CST_Theta_AP_mean SLF_Theta_AP_mean ; CC_Theta_AP_stdev CR_Theta_AP_stdev CST_Theta_AP_stdev SLF_Theta_AP_stdev]
AP_SPERC = [CC_Sperc_AP CR_Sperc_AP CST_Sperc_AP SLF_Sperc_AP]
AP_FPERC = [CC_Fperc_AP CR_Fperc_AP CST_Fperc_AP SLF_Fperc_AP]

LR_VSS = [CC_VSS_LR_mean CR_VSS_LR_mean CST_VSS_LR_mean SLF_VSS_LR_mean ; CC_VSS_LR_stdev CR_VSS_LR_stdev CST_VSS_LR_stdev SLF_VSS_LR_stdev]
LR_THETA = [CC_Theta_LR_mean CR_Theta_LR_mean CST_Theta_LR_mean SLF_Theta_LR_mean ; CC_Theta_LR_stdev CR_Theta_LR_stdev CST_Theta_LR_stdev SLF_Theta_LR_stdev]
LR_SPERC = [CC_Sperc_LR CR_Sperc_LR CST_Sperc_LR SLF_Sperc_LR]
LR_FPERC = [CC_Fperc_LR CR_Fperc_LR CST_Fperc_LR SLF_Fperc_LR]

SI_VSS = [CC_VSS_SI_mean CR_VSS_SI_mean CST_VSS_SI_mean SLF_VSS_SI_mean ; CC_VSS_SI_stdev CR_VSS_SI_stdev CST_VSS_SI_stdev SLF_VSS_SI_stdev]
SI_THETA = [CC_Theta_SI_mean CR_Theta_SI_mean CST_Theta_SI_mean SLF_Theta_SI_mean ; CC_Theta_SI_stdev CR_Theta_SI_stdev CST_Theta_SI_stdev SLF_Theta_SI_stdev]
SI_SPERC = [CC_Sperc_SI CR_Sperc_SI CST_Sperc_SI SLF_Sperc_SI]
SI_FPERC = [CC_Fperc_SI CR_Fperc_SI CST_Fperc_SI SLF_Fperc_SI]

