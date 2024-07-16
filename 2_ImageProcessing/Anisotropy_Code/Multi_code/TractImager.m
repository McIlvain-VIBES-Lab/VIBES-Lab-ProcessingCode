load('dtiall.mat');
%% Load in the Tracts
cd Tracts
CC_ = load_nii('tractsCC.nii'); CC = repmat(CC_.img>.6,[1 1 1 3]); CCT = CC_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CR_ = load_nii('tractsCR.nii'); CR = repmat(CR_.img>.6,[1 1 1 3]); CRT = CR_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
CST_ = load_nii('tractsCCS.nii'); CST = repmat(CST_.img>.6, [1 1 1 3]); CSTT = CST_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
SLF_ = load_nii('tractsSLF.nii'); SLF = repmat(SLF_.img>.6, [1 1 1 3]); SLFT = SLF_.img.*(abs(dtiall(:,:,:,1))>0)>.6;
cd ..

Tracts = CCT+CRT+CSTT+SLFT;
%% Load in AP F/S percentage 
% Means and Standard Devs
cd MRE; cd AP;
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd .. 
Ussum = sum(abs(U_si),5); Ufsum = sum(abs(U_fi),5);
Uf = sum(Ufsum,4); Us = sum(Ussum,4); U = Us+Uf;
Usperc_AP = Us./U; Ufperc_AP = Uf./U;

Us40AP = (Usperc_AP>.40).*Tracts; Us50AP = (Usperc_AP>.5).*Tracts; Us60AP = (Usperc_AP>.60).*Tracts;
Uf40AP = (Ufperc_AP>.40).*Tracts; Uf50AP = (Ufperc_AP>.5).*Tracts; Uf60AP = (Ufperc_AP>.60).*Tracts;

%% Load in LR F/S percentage 
% Means and Standard Devs
cd LR;
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd .. 

Ussum = sum(abs(U_si),5); Ufsum = sum(abs(U_fi),5);
Uf = sum(Ufsum,4); Us = sum(Ussum,4); U = Us+Uf;
Usperc_LR = Us./U; Ufperc_LR = Uf./U;

Us40LR = (Usperc_LR>.40).*Tracts; Us50LR = (Usperc_LR>.5).*Tracts; Us60LR = (Usperc_LR>.60).*Tracts;
Uf40LR = (Ufperc_LR>.40).*Tracts; Uf50LR = (Ufperc_LR>.5).*Tracts; Uf60LR = (Ufperc_LR>.60).*Tracts;

%% Load in SI F/S percentage 
% Means and Standard Devs
cd SI;
cd Iterative_Outputs; load Outputs.mat; load Waveprop_BMR.mat; cd .. % S_vox_U F_vox_U Us_norm Uf_norm 
cd .. 

Ussum = sum(abs(U_si),5); Ufsum = sum(abs(U_fi),5);
Uf = sum(Ufsum,4); Us = sum(Ussum,4); U = Us+Uf;
Usperc_LR = Us./U; Ufperc_LR = Uf./U;
 
Us40SI = (Usperc_SI>.40).*Tracts; Us50SI = (Usperc_SI>.5).*Tracts; Us60SI = (Usperc_SI>.6).*Tracts;
Uf40SI = (Ufperc_SI>.40).*Tracts; Uf50SI = (Ufperc_SI>.5).*Tracts; Uf60SI = (Ufperc_SI>.6).*Tracts;

%% Figure Creation

figure;im(Us40AP); figure;im(Us50AP); figure;im(Us60AP); figure;im(Uf40AP); figure;im(Uf50AP); figure;im(Uf60AP);
figure;im(Us40LR); figure;im(Us50LR); figure;im(Us60LR); figure;im(Uf40LR); figure;im(Uf50LR); figure;im(Uf60LR);
figure;im(Us40SI); figure;im(Us50SI); figure;im(Us60SI); figure;im(Uf40SI); figure;im(Uf50SI); figure;im(Uf60SI);

cd ..
