load('dtiall.mat');
cd MRE
cd Tracts
CC = load_nii('tractsCC.nii'); CC = CC.img>.6;
CR = load_nii('tractsCR.nii'); CR = CR.img>.6;
CCS = load_nii('tractsCCS.nii'); CCS = CCS.img>.6;
SLF = load_nii('tractsSLF.nii'); SLF = SLF.img>.6;
cd ..

cd AP
load Waveprop_BMR.mat
cd PostInv
load Mu.mat
AP_CC_Mu = Mu.*CC.*dtiall(:,:,:,1)>0; AP_CR_Mu = Mu.*CR.*dtiall(:,:,:,1)>0; AP_CCS_Mu = Mu.*CCS.*dtiall(:,:,:,1)>0; AP_SLF_Mu = Mu.*SLF.*dtiall(:,:,:,1)>0;
AP_CC_theta = theta2.*CC; AP_CR_theta = theta2.*CR; AP_CCS_theta = theta2.*CCS; AP_SLF_theta = theta2.*SLF;
AP_CC_Us = abs(Usd_bar).*CC; AP_CR_Us = abs(Usd_bar).*CR; AP_CCS_Us = abs(Usd_bar).*CCS; AP_SLF_Us = abs(Usd_bar).*SLF;
AP_CC_Uf = abs(Ufd_bar).*CC; AP_CR_Uf = abs(Ufd_bar).*CR; AP_CCS_Uf = abs(Ufd_bar).*CCS; AP_SLF_Uf = abs(Ufd_bar).*SLF;

AP_CC_Mu(isnan(AP_CC_Mu)) = 0; AP_CR_Mu(isnan(AP_CR_Mu)) = 0; AP_CCS_Mu(isnan(AP_CCS_Mu)) = 0; AP_SLF_Mu(isnan(AP_SLF_Mu)) = 0;
AP_CC_theta(isnan(AP_CC_theta)) = 0; AP_CR_Mu(isnan(AP_CR_theta)) = 0; AP_CCS_theta(isnan(AP_CCS_theta)) = 0; AP_SLF_theta(isnan(AP_SLF_theta)) = 0;
AP_CC_Us(isnan(AP_CC_Us)) = 0; AP_CR_Us(isnan(AP_CR_Us)) = 0; AP_CCS_Us(isnan(AP_CCS_Us)) = 0; AP_SLF_Us(isnan(AP_SLF_Us)) = 0; 
AP_CC_Uf(isnan(AP_CC_Uf)) = 0; AP_CR_Uf(isnan(AP_CR_Uf)) = 0; AP_CCS_Uf(isnan(AP_CCS_Uf)) = 0; AP_SLF_Uf(isnan(AP_SLF_Uf)) = 0; 

AP_CC_Mu0 = AP_CC_Mu(AP_CC_Mu>0);AP_CR_Mu0 = AP_CR_Mu(AP_CR_Mu>0); AP_CCS_Mu0 = AP_CCS_Mu(AP_CCS_Mu>0); AP_SLF_Mu0 = AP_SLF_Mu(AP_SLF_Mu>0);
AP_CC_theta0 = AP_CC_theta(AP_CC_theta>0);AP_CR_theta0 = AP_CR_theta(AP_CR_theta>0); AP_CCS_theta0 = AP_CCS_theta(AP_CCS_theta>0); AP_SLF_theta0 = AP_SLF_theta(AP_SLF_theta>0);
AP_CC_Us0 = AP_CC_Us(AP_CC_Us>0); AP_CR_Us0 = AP_CR_Us(AP_CR_Us>0); AP_CCS_Us0 = AP_CCS_Us(AP_CCS_Us>0); AP_SLF_Us0 = AP_SLF_Us(AP_SLF_Us>0);
AP_CC_Uf0 = AP_CC_Uf(AP_CC_Uf>0); AP_CR_Uf0 = AP_CR_Uf(AP_CR_Uf>0); AP_CCS_Uf0 = AP_CCS_Uf(AP_CCS_Uf>0); AP_SLF_Uf0 = AP_SLF_Uf(AP_SLF_Uf>0);

AP_CC_Mu_mean = mean(AP_CC_Mu0); AP_CR_Mu_mean = mean(AP_CR_Mu0); AP_CCS_Mu_mean = mean(AP_CCS_Mu0); AP_SLF_Mu_mean = mean(AP_SLF_Mu0);   
AP_CC_theta_mean = mean(AP_CC_theta0); AP_CR_theta_mean = mean(AP_CR_theta0); AP_CCS_theta_mean = mean(AP_CCS_theta0); AP_SLF_theta_mean = mean(AP_SLF_theta0);  
AP_CC_Us_mean = mean(AP_CC_Us0); AP_CR_Us_mean = mean(AP_CR_Us0); AP_CCS_Us_mean = mean(AP_CCS_Us0); AP_SLF_Us_mean = mean(AP_SLF_Us0);   
AP_CC_Uf_mean = mean(AP_CC_Uf0); AP_CR_Uf_mean = mean(AP_CR_Uf0); AP_CCS_Uf_mean = mean(AP_CCS_Uf0); AP_SLF_Uf_mean = mean(AP_SLF_Uf0);   
cd ..
cd ..

cd LR
load Waveprop_BMR.mat
cd PostInv
load Mu.mat
LR_CC_Mu = Mu.*CC; LR_CR_Mu = Mu.*CR; LR_CCS_Mu = Mu.*CCS; LR_SLF_Mu = Mu.*SLF;
LR_CC_theta = theta2.*CC; LR_CR_theta = theta2.*CR; LR_CCS_theta = theta2.*CCS; LR_SLF_theta = theta2.*SLF;
LR_CC_Us = abs(Usd_bar).*CC; LR_CR_Us = abs(Usd_bar).*CR; LR_CCS_Us = abs(Usd_bar).*CCS; LR_SLF_Us = abs(Usd_bar).*SLF;
LR_CC_Uf = abs(Ufd_bar).*CC; LR_CR_Uf = abs(Ufd_bar).*CR; LR_CCS_Uf = abs(Ufd_bar).*CCS; LR_SLF_Uf = abs(Ufd_bar).*SLF;

LR_CC_Mu(isnan(LR_CC_Mu)) = 0; LR_CR_Mu(isnan(LR_CR_Mu)) = 0; LR_CCS_Mu(isnan(LR_CCS_Mu)) = 0; LR_SLF_Mu(isnan(LR_SLF_Mu)) = 0;
LR_CC_theta(isnan(LR_CC_theta)) = 0; LR_CR_Mu(isnan(LR_CR_theta)) = 0; LR_CCS_theta(isnan(LR_CCS_theta)) = 0; LR_SLF_theta(isnan(LR_SLF_theta)) = 0;
LR_CC_Us(isnan(LR_CC_Us)) = 0; LR_CR_Us(isnan(LR_CR_Us)) = 0; LR_CCS_Us(isnan(LR_CCS_Us)) = 0; LR_SLF_Us(isnan(LR_SLF_Us)) = 0; 
LR_CC_Uf(isnan(LR_CC_Uf)) = 0; LR_CR_Uf(isnan(LR_CR_Uf)) = 0; LR_CCS_Uf(isnan(LR_CCS_Uf)) = 0; LR_SLF_Uf(isnan(LR_SLF_Uf)) = 0; 

LR_CC_Mu0 = LR_CC_Mu(LR_CC_Mu>0); LR_CR_Mu0 = LR_CR_Mu(LR_CR_Mu>0); LR_CCS_Mu0 = LR_CCS_Mu(LR_CCS_Mu>0); LR_SLF_Mu0 = LR_SLF_Mu(LR_SLF_Mu>0);
LR_CC_theta0 = LR_CC_theta(LR_CC_theta>0); LR_CR_theta0 =LR_CR_theta(LR_CR_theta>0); LR_CCS_theta0 = LR_CCS_theta(LR_CCS_theta>0); LR_SLF_theta0 = LR_SLF_theta(LR_SLF_theta>0);
LR_CC_Us0 = LR_CC_Us(LR_CC_Us>0); LR_CR_Us0 = LR_CR_Us(LR_CR_Us>0); LR_CCS_Us0 = LR_CCS_Us(LR_CCS_Us>0); LR_SLF_Us0 = LR_SLF_Us(LR_SLF_Us>0);
LR_CC_Uf0 = LR_CC_Uf(LR_CC_Uf>0); LR_CR_Uf0 = LR_CR_Uf(LR_CR_Uf>0); LR_CCS_Uf0 = LR_CCS_Uf(LR_CCS_Uf>0); LR_SLF_Uf0 = LR_SLF_Uf(LR_SLF_Uf>0);

LR_CC_Mu_mean = mean(LR_CC_Mu0); LR_CR_Mu_mean = mean(LR_CR_Mu0); LR_CCS_Mu_mean = mean(LR_CCS_Mu0); LR_SLF_Mu_mean = mean(LR_SLF_Mu0);   
LR_CC_theta_mean = mean(LR_CC_theta0); LR_CR_theta_mean = mean(LR_CR_theta0); LR_CCS_theta_mean = mean(LR_CCS_theta0); LR_SLF_theta_mean = mean(LR_SLF_theta0);  
LR_CC_Us_mean = mean(LR_CC_Us0); LR_CR_Us_mean = mean(LR_CR_Us0); LR_CCS_Us_mean = mean(LR_CCS_Us0); LR_SLF_Us_mean = mean(LR_SLF_Us0);   
LR_CC_Uf_mean = mean(LR_CC_Uf0); LR_CR_Uf_mean = mean(LR_CR_Uf0); LR_CCS_Uf_mean = mean(LR_CCS_Uf0); LR_SLF_Uf_mean = mean(LR_SLF_Uf0);   

cd ..
cd ..

cd SI
load Waveprop_BMR.mat
cd PostInv
load Mu.mat
SI_CC_Mu = Mu.*CC; SI_CR_Mu = Mu.*CR; SI_CCS_Mu = Mu.*CCS; SI_SLF_Mu = Mu.*SLF;
SI_CC_theta = theta2.*CC; SI_CR_theta = theta2.*CR; SI_CCS_theta = theta2.*CCS; SI_SLF_theta = theta2.*SLF;
SI_CC_Us = abs(Usd_bar).*CC; SI_CR_Us = abs(Usd_bar).*CR; SI_CCS_Us = abs(Usd_bar).*CCS; SI_SLF_Us = abs(Usd_bar).*SLF;
SI_CC_Uf = abs(Ufd_bar).*CC; SI_CR_Uf = abs(Ufd_bar).*CR; SI_CCS_Uf = abs(Ufd_bar).*CCS; SI_SLF_Uf = abs(Ufd_bar).*SLF;

SI_CC_Mu(isnan(SI_CC_Mu)) = 0; SI_CR_Mu(isnan(SI_CR_Mu)) = 0; SI_CCS_Mu(isnan(SI_CCS_Mu)) = 0; SI_SLF_Mu(isnan(SI_SLF_Mu)) = 0;
SI_CC_theta(isnan(SI_CC_theta)) = 0; SI_CR_Mu(isnan(SI_CR_theta)) = 0; SI_CCS_theta(isnan(SI_CCS_theta)) = 0; SI_SLF_theta(isnan(SI_SLF_theta)) = 0;
SI_CC_Us(isnan(SI_CC_Us)) = 0; SI_CR_Us(isnan(SI_CR_Us)) = 0; SI_CCS_Us(isnan(SI_CCS_Us)) = 0; SI_SLF_Us(isnan(SI_SLF_Us)) = 0; 
SI_CC_Uf(isnan(SI_CC_Uf)) = 0; SI_CR_Uf(isnan(SI_CR_Uf)) = 0; SI_CCS_Uf(isnan(SI_CCS_Uf)) = 0; SI_SLF_Uf(isnan(SI_SLF_Uf)) = 0; 

SI_CC_Mu0 = SI_CC_Mu(SI_CC_Mu>0); SI_CR_Mu0 = SI_CR_Mu(SI_CR_Mu>0); SI_CCS_Mu0 = SI_CCS_Mu(SI_CCS_Mu>0); SI_SLF_Mu0 = SI_SLF_Mu(SI_SLF_Mu>0);
SI_CC_theta0 = SI_CC_theta(SI_CC_theta>0); SI_CR_theta0 = SI_CR_theta(SI_CR_theta>0); SI_CCS_theta0 = SI_CCS_theta(SI_CCS_theta>0); SI_SLF_theta0 = SI_SLF_theta(SI_SLF_theta>0);
SI_CC_Us0 = SI_CC_Us(SI_CC_Us>0); SI_CR_Us0 = SI_CR_Us(SI_CR_Us>0); SI_CCS_Us0 = SI_CCS_Us(SI_CCS_Us>0); SI_SLF_Us0 = SI_SLF_Us(SI_SLF_Us>0);
SI_CC_Uf0 = SI_CC_Uf(SI_CC_Uf>0); SI_CR_Uf0 = SI_CR_Uf(SI_CR_Uf>0); SI_CCS_Uf0 = SI_CCS_Uf(SI_CCS_Uf>0); SI_SLF_Uf0 = SI_SLF_Uf(SI_SLF_Uf>0);

SI_CC_Mu_mean = mean(SI_CC_Mu0); SI_CR_Mu_mean = mean(SI_CR_Mu0); SI_CCS_Mu_mean = mean(SI_CCS_Mu0); SI_SLF_Mu_mean = mean(SI_SLF_Mu0);   
SI_CC_theta_mean = mean(SI_CC_theta0); SI_CR_theta_mean = mean(SI_CR_theta0); SI_CCS_theta_mean = mean(SI_CCS_theta0); SI_SLF_theta_mean = mean(SI_SLF_theta0);  
SI_CC_Us_mean = mean(SI_CC_Us0); SI_CR_Us_mean = mean(SI_CR_Us0); SI_CCS_Us_mean = mean(SI_CCS_Us0); SI_SLF_Us_mean = mean(SI_SLF_Us0);   
SI_CC_Uf_mean = mean(SI_CC_Uf0); SI_CR_Uf_mean = mean(SI_CR_Uf0); SI_CCS_Uf_mean = mean(SI_CCS_Uf0); SI_SLF_Uf_mean = mean(SI_SLF_Uf0);   
cd ..
cd ..
cd ..

AP_CC = [AP_CC_Mu_mean, AP_CC_theta_mean, AP_CC_Us_mean, AP_CC_Uf_mean; std2(AP_CC_Mu0), std2(AP_CC_theta0), std2(AP_CC_Us0), std2(AP_CC_Uf0)] 
AP_CR = [AP_CR_Mu_mean, AP_CR_theta_mean, AP_CR_Us_mean, AP_CR_Uf_mean; std2(AP_CR_Mu0), std2(AP_CR_theta0), std2(AP_CR_Us0), std2(AP_CR_Uf0)]
AP_CCS = [AP_CCS_Mu_mean, AP_CCS_theta_mean, AP_CCS_Us_mean, AP_CCS_Uf_mean; std2(AP_CCS_Mu0), std2(AP_CCS_theta0), std2(AP_CCS_Us0), std2(AP_CCS_Uf0)] 
AP_SLF = [AP_SLF_Mu_mean, AP_SLF_theta_mean, AP_SLF_Us_mean, AP_SLF_Uf_mean; std2(AP_SLF_Mu0), std2(AP_SLF_theta0), std2(AP_SLF_Us0), std2(AP_SLF_Uf0)] 
LR_CC = [LR_CC_Mu_mean, LR_CC_theta_mean, LR_CC_Us_mean, LR_CC_Uf_mean; std2(LR_CC_Mu0), std2(LR_CC_theta0), std2(LR_CC_Us0), std2(LR_CC_Uf0)]
LR_CR = [LR_CR_Mu_mean, LR_CR_theta_mean, LR_CR_Us_mean, LR_CR_Uf_mean; std2(LR_CR_Mu0), std2(LR_CR_theta0), std2(LR_CR_Us0), std2(LR_CR_Uf0)]
LR_CCS = [LR_CCS_Mu_mean, LR_CCS_theta_mean, LR_CCS_Us_mean, LR_CCS_Uf_mean; std2(LR_CCS_Mu0), std2(LR_CCS_theta0), std2(LR_CCS_Us0), std2(LR_CCS_Uf0)] 
LR_SLF = [LR_SLF_Mu_mean, LR_SLF_theta_mean, LR_SLF_Us_mean, LR_SLF_Uf_mean; std2(LR_SLF_Mu0), std2(LR_SLF_theta0), std2(LR_SLF_Us0), std2(LR_SLF_Uf0)] 
SI_CC = [SI_CC_Mu_mean, SI_CC_theta_mean, SI_CC_Us_mean, SI_CC_Uf_mean; std(SI_CC_Mu0), std(SI_CC_theta0), std(SI_CC_Us0), std(SI_CC_Uf0)]
SI_CR = [SI_CR_Mu_mean, SI_CR_theta_mean, SI_CR_Us_mean, SI_CR_Uf_mean; std(SI_CR_Mu0), std(SI_CR_theta0), std(SI_CR_Us0), std(SI_CR_Uf0)]
SI_CCS = [SI_CCS_Mu_mean, SI_CCS_theta_mean, SI_CCS_Us_mean, SI_CCS_Uf_mean; std(SI_CCS_Mu0), std(SI_CCS_theta0), std(SI_CCS_Us0), std(SI_CCS_Uf0)] 
SI_SLF = [SI_SLF_Mu_mean, SI_SLF_theta_mean, SI_SLF_Us_mean, SI_SLF_Uf_mean; std(SI_SLF_Mu0), std(SI_SLF_theta0), std(SI_SLF_Us0), std(SI_SLF_Uf0)] 
voxels = [size(AP_CC_Mu0), size(AP_CR_Mu0), size(AP_SLF_Mu0), size(AP_CCS_Mu0)]