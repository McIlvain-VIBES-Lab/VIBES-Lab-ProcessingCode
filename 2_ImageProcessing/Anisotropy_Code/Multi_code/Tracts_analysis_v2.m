load('dtiall.mat');

cd Tracts
CC_ = load_nii('tractsCC.nii'); CC = repmat(CC_.img>.6,[1 1 1 3]); CCT = CC_.img>.6;
CR_ = load_nii('tractsCR.nii'); CR = repmat(CR_.img>.6,[1 1 1 3]); CRT = CR_.img>.6;
CCS_ = load_nii('tractsCCS.nii'); CCS = repmat(CCS_.img>.6, [1 1 1 3]); CCST = CCS_.img>.6;
SLF_ = load_nii('tractsSLF.nii'); SLF = repmat(SLF_.img>.6, [1 1 1 3]); SLFT = SLF_.img>.6;
cd ..

cd MRE
cd AP
cd Iterative_Outputs
load Waveprop_BMR.mat
load Outputs.mat
Fvoxel = repmat(F_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Fvox = F_vox.*(dtiall(:,:,:,1))>0;
Svoxel = repmat(S_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Svox = S_vox.*(dtiall(:,:,:,1))>0;
Fast = Nf_norm.*Fvoxel; Slow = Ns_norm.*Svoxel;
CCAP_FA = CC.*Fvoxel.*dtiall; CRAP_FA = CR.*Fvoxel.*dtiall; CCSAP_FA = CCS.*Fvoxel.*dtiall; SLFAP_FA = SLF.*Fvoxel.*dtiall;
CCAP_SA = CC.*Svoxel.*dtiall; CRAP_SA = CR.*Svoxel.*dtiall; CCSAP_SA = CCS.*Svoxel.*dtiall; SLFAP_SA = SLF.*Svoxel.*dtiall;
CCAP_FN = CC.*Fast; CRAP_FN = CR.*Fast; CCSAP_FN = CCS.*Fast; SLFAP_FN = SLF.*Fast;
CCAP_SN = CC.*Slow; CRAP_SN = CR.*Slow; CCSAP_SN = CCS.*Slow; SLFAP_SN = SLF.*Slow;
CCAP_FTheta = CCT.*Theta_fast.*Fvox; CRAP_FTheta = CRT.*Theta_fast.*Fvox; 
CCSAP_FTheta = CCST.*Theta_fast.*Fvox; SLFAP_FTheta = SLFT.*Theta_fast.*Fvox;
CCAP_STheta = CCT.*Theta_slow.*Svox; CRAP_STheta = CRT.*Theta_slow.*Svox; 
CCSAP_STheta = CCST.*Theta_slow.*Svox; SLFAP_STheta = SLFT.*Theta_slow.*Svox;
cd ..
cd ..

cd LR
cd Iterative_Outputs
load Waveprop_BMR.mat
load Outputs.mat
Fvoxel = repmat(F_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Fvox = F_vox.*(dtiall(:,:,:,1))>0;
Svoxel = repmat(S_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Svox = S_vox.*(dtiall(:,:,:,1))>0;
Fast = Nf_norm.*Fvoxel; Slow = Ns_norm.*Svoxel;
CCLR_FA = CC.*Fvoxel.*dtiall; CRLR_FA = CR.*Fvoxel.*dtiall; CCSLR_FA = CCS.*Fvoxel.*dtiall; SLFLR_FA = SLF.*Fvoxel.*dtiall;
CCLR_SA = CC.*Svoxel.*dtiall; CRLR_SA = CR.*Svoxel.*dtiall; CCSLR_SA = CCS.*Svoxel.*dtiall; SLFLR_SA = SLF.*Svoxel.*dtiall;
CCLR_FN = CC.*Fast; CRLR_FN = CR.*Fast; CCSLR_FN = CCS.*Fast; SLFLR_FN = SLF.*Fast;
CCLR_SN = CC.*Slow; CRLR_SN = CR.*Slow; CCSLR_SN = CCS.*Slow; SLFLR_SN = SLF.*Slow;
CCLR_FTheta = CCT.*Theta_fast.*Fvox; CRLR_FTheta = CRT.*Theta_fast.*Fvox; 
CCSLR_FTheta = CCST.*Theta_fast.*Fvox; SLFLR_FTheta = SLFT.*Theta_fast.*Fvox;
CCLR_STheta = CCT.*Theta_slow.*Svox; CRLR_STheta = CRT.*Theta_slow.*Svox; 
CCSLR_STheta = CCST.*Theta_slow.*Svox; SLFLR_STheta = SLFT.*Theta_slow.*Svox;
cd ..
cd ..

cd SI
cd Iterative_Outputs
load Waveprop_BMR.mat
load Outputs.mat
Fvoxel = repmat(F_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Fvox = F_vox.*(dtiall(:,:,:,1))>0;
Svoxel = repmat(S_vox.*(dtiall(:,:,:,1))>0,[1 1 1 3]);
Svox = S_vox.*(dtiall(:,:,:,1))>0;
Fast = Nf_norm.*Fvoxel; Slow = Ns_norm.*Svoxel;
CCSI_FA = CC.*Fvoxel.*dtiall; CRSI_FA = CR.*Fvoxel.*dtiall; CCSSI_FA = CCS.*Fvoxel.*dtiall; SLFSI_FA = SLF.*Fvoxel.*dtiall;
CCSI_SA = CC.*Svoxel.*dtiall; CRSI_SA = CR.*Svoxel.*dtiall; CCSSI_SA = CCS.*Svoxel.*dtiall; SLFSI_SA = SLF.*Svoxel.*dtiall;
CCSI_FN = CC.*Fast; CRSI_FN = CR.*Fast; CCSSI_FN = CCS.*Fast; SLFSI_FN = SLF.*Fast;
CCSI_SN = CC.*Slow; CRSI_SN = CR.*Slow; CCSSI_SN = CCS.*Slow; SLFSI_SN = SLF.*Slow;
CCSI_FTheta = CCT.*Theta_fast.*Fvox; CRSI_FTheta = CRT.*Theta_fast.*Fvox; 
CCSSI_FTheta = CCST.*Theta_fast.*Fvox; SLFSI_FTheta = SLFT.*Theta_fast.*Fvox;
CCSI_STheta = CCT.*Theta_slow.*Svox; CRSI_STheta = CRT.*Theta_slow.*Svox;
CCSSI_STheta = CCST.*Theta_slow.*Svox; SLFSI_STheta = SLFT.*Theta_slow.*Svox;
cd ..
cd .. 
cd ..

APFast_A = [mean(abs(CCAP_FA(CCAP_FA>0))); mean(abs(CRAP_FA(CRAP_FA>0))); mean(abs(CCSAP_FA(CCSAP_FA>0))); mean(abs(SLFAP_FA(SLFAP_FA>0)))]
APSlow_A = [mean(abs(CCAP_SA(CCAP_SA>0))); mean(abs(CRAP_SA(CRAP_SA>0))); mean(abs(CCSAP_SA(CCSAP_SA>0))); mean(abs(SLFAP_SA(SLFAP_SA>0)))]
APFast_N = [mean(abs(CCAP_FN(CCAP_FN>0))); mean(abs(CRAP_FN(CRAP_FN>0))); mean(abs(CCSAP_FN(CCSAP_FN>0))); mean(abs(SLFAP_FN(SLFAP_FN>0)))]
APSlow_N = [mean(abs(CCAP_SN(CCAP_SN>0))); mean(abs(CRAP_SN(CRAP_SN>0))); mean(abs(CCSAP_SN(CCSAP_SN>0))); mean(abs(SLFAP_SN(SLFAP_SN>0)))]
APFast_Theta = [mean(abs(CCAP_FTheta(CCAP_FTheta>0))); mean(abs(CRAP_FTheta(CRAP_FTheta>0))); mean(abs(CCSAP_FTheta(CCSAP_FTheta>0))); mean(abs(SLFAP_FTheta(SLFAP_FTheta>0)))]
APSlow_Theta = [mean(abs(CCAP_STheta(CCAP_STheta>0))); mean(abs(CRAP_STheta(CRAP_STheta>0))); mean(abs(CCSAP_STheta(CCSAP_STheta>0))); mean(abs(SLFAP_STheta(SLFAP_STheta>0)))]
LRFast_A = [mean(abs(CCLR_FA(CCLR_FA>0))); mean(abs(CRLR_FA(CRLR_FA>0))); mean(abs(CCSLR_FA(CCSLR_FA>0))); mean(abs(SLFLR_FA(SLFLR_FA>0)))]
LRSlow_A = [mean(abs(CCLR_SA(CCLR_SA>0))); mean(abs(CRLR_SA(CRLR_SA>0))); mean(abs(CCSLR_SA(CCSLR_SA>0))); mean(abs(SLFLR_SA(SLFLR_SA>0)))]
LRFast_N = [mean(abs(CCLR_FN(CCLR_FN>0))); mean(abs(CRLR_FN(CRLR_FN>0))); mean(abs(CCSLR_FN(CCSLR_FN>0))); mean(abs(SLFLR_FN(SLFLR_FN>0)))]
LRSlow_N = [mean(abs(CCLR_SN(CCLR_SN>0))); mean(abs(CRLR_SN(CRLR_SN>0))); mean(abs(CCSLR_SN(CCSLR_SN>0))); mean(abs(SLFLR_SN(SLFLR_SN>0)))]
LRFast_Theta = [mean(abs(CCLR_FTheta(CCLR_FTheta>0))); mean(abs(CRLR_FTheta(CRLR_FTheta>0))); mean(abs(CCSLR_FTheta(CCSLR_FTheta>0))); mean(abs(SLFLR_FTheta(SLFLR_FTheta>0)))]
LRSlow_Theta = [mean(abs(CCLR_STheta(CCLR_STheta>0))); mean(abs(CRLR_STheta(CRLR_STheta>0))); mean(abs(CCSLR_STheta(CCSLR_STheta>0))); mean(abs(SLFLR_STheta(SLFLR_STheta>0)))]
SIFast_A = [mean(abs(CCSI_FA(CCSI_FA>0))); mean(abs(CRSI_FA(CRSI_FA>0))); mean(abs(CCSSI_FA(CCSSI_FA>0))); mean(abs(SLFSI_FA(SLFSI_FA>0)))]
SISlow_A = [mean(abs(CCSI_SA(CCSI_SA>0))); mean(abs(CRSI_SA(CRSI_SA>0))); mean(abs(CCSSI_SA(CCSSI_SA>0))); mean(abs(SLFSI_SA(SLFSI_SA>0)))]
SIFast_N = [mean(abs(CCSI_FN(CCSI_FN>0))); mean(abs(CRSI_FN(CRSI_FN>0))); mean(abs(CCSSI_FN(CCSSI_FN>0))); mean(abs(SLFSI_FN(SLFSI_FN>0)))]
SISlow_N = [mean(abs(CCSI_SN(CCSI_SN>0))); mean(abs(CRSI_SN(CRSI_SN>0))); mean(abs(CCSSI_SN(CCSSI_SN>0))); mean(abs(SLFSI_SN(SLFSI_SN>0)))]
SIFast_Theta = [mean(abs(CCSI_FTheta(CCSI_FTheta>0))); mean(abs(CRSI_FTheta(CRSI_FTheta>0))); mean(abs(CCSSI_FTheta(CCSSI_FTheta>0))); mean(abs(SLFSI_FTheta(SLFSI_FTheta>0)))]
SISlow_Theta = [mean(abs(CCSI_STheta(CCSI_STheta>0))); mean(abs(CRSI_STheta(CRSI_STheta>0))); mean(abs(CCSSI_STheta(CCSSI_STheta>0))); mean(abs(SLFSI_STheta(SLFSI_STheta>0)))]
