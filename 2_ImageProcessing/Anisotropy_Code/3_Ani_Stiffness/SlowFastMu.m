cd /Volumes/CLJ-001-1/Group/drsmitty/ScanData/Multi/170222
cd Multiexcite_01; cd MRE; cd AP
load Quant_MUFS_AP.mat
cd ..; cd LR
load Quant_MUFS_LR.mat
cd ..


FAmin = .25;
FAmax = .85;

R_fa = zeros(13,43);
R_fa_AP = zeros(13,61);
R_fa_LR = zeros(13,61);

for Tract = 1
    Tract
    [Rfa,Rfa_AP,Rfa_LR,FAall] = FA_correlation(FAmin,FAmax,Tract,FAtractlistAP,FAtractlistLR,MutotlistAP,UperctotAP,MutotlistLR,UperctotLR,AmplistAP,AmplistLR);
    R_fa(Tract,1:length(Rfa)) = abs(Rfa);
    R_fa_AP(Tract,1:length(Rfa)) = abs(Rfa_AP);
    R_fa_LR(Tract,1:length(Rfa)) = abs(Rfa_LR);
    
    R(1,Tract) = max(abs(Rfa));
    R(2,Tract) = max(abs(Rfa_AP));
    R(3,Tract) = max(abs(Rfa_LR));    
    
    Rfaall = FAall(Rfa == max(abs(Rfa)));
    RfaLR = FAall(Rfa_LR == max(abs(Rfa_LR)));
    RfaAP = FAall(Rfa_AP == max(abs(Rfa_AP))); 
    
    FAR(1,Tract) = Rfaall(1);
    FAR(2,Tract) = RfaAP(1);
    FAR(3,Tract) = RfaLR(1);    
end

RCC(1,:) = R_fa(1,:); RCC(2,:) = R_fa(6,:); RCC(3,:) = R_fa(13,:);
RCC_AP(1,:) = R_fa_AP(1,:); RCC_AP(2,:) = R_fa_AP(6,:); RCC_AP(3,:) = R_fa_AP(13,:);
RCC_LR(1,:) = R_fa_LR(1,:); RCC_LR(2,:) = R_fa_LR(6,:); RCC_LR(3,:) = R_fa_LR(13,:);
figure;plot(FAall,RCC)
figure;plot(FAall,RCC_AP)
figure;plot(FAall,RCC_LR)

%[Slope1,inter1,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval1] = TractSF(FAR(1,:),MutotlistAP,MutotlistLR,UperctotAP,UperctotLR,FAtractlistAP,FAtractlistLR,AmplistAP,AmplistLR);

FAtractlistAP1 = FAtractlistAP;
FAtractlistLR1 = FAtractlistLR;
UperctotAP1 = UperctotAP;
UperctotLR1 = UperctotLR;
MutotlistAP1 = MutotlistAP;
MutotlistLR1 = MutotlistLR;
AmplistAP1 = AmplistAP;
AmplistLR1 = AmplistLR;

cd /Volumes/CLJ-001-1/Group/drsmitty/ScanData/Multi/170222

cd Multiexcite_02; cd MRE; cd AP
load Quant_MUFS_AP.mat
cd ..; cd LR
load Quant_MUFS_LR.mat
cd ..
FAmin = .25;
FAmax = .85;
R_fa = zeros(13,43);
R_fa_AP = zeros(13,61);
R_fa_LR = zeros(13,61);

for Tract = 1:13
    Tract
    [Rfa,Rfa_AP,Rfa_LR,FAall] = FA_correlation(FAmin,FAmax,Tract,FAtractlistAP,FAtractlistLR,MutotlistAP,UperctotAP,MutotlistLR,UperctotLR,AmplistAP,AmplistLR);
    R_fa(Tract,1:length(Rfa)) = abs(Rfa);
    R_fa_AP(Tract,1:length(Rfa)) = abs(Rfa_AP);
    R_fa_LR(Tract,1:length(Rfa)) = abs(Rfa_LR);
    
    R(1,Tract) = max(abs(Rfa));
    R(2,Tract) = max(abs(Rfa_AP));
    R(3,Tract) = max(abs(Rfa_LR));    
    
    Rfaall = FAall(Rfa == max(abs(Rfa)));
    RfaLR = FAall(Rfa_LR == max(abs(Rfa_LR)));
    RfaAP = FAall(Rfa_AP == max(abs(Rfa_AP)));
    
    FAR(1,Tract) = Rfaall(1);
    FAR(2,Tract) = RfaAP(1);
    FAR(3,Tract) = RfaLR(1);    
end

RCC(1,:) = R_fa(1,:); RCC(2,:) = R_fa(6,:); RCC(3,:) = R_fa(13,:);
RCC_AP(1,:) = R_fa_AP(1,:); RCC(2,:) = R_fa_AP(6,:); RCC(3,:) = R_fa_AP(13,:);
RCC_LR(1,:) = R_fa_LR(1,:); RCC(2,:) = R_fa_LR(6,:); RCC(3,:) = R_fa_LR(13,:);
figure;plot(FAall,RCC)
figure;plot(FAall,RCC_AP)
figure;plot(FAall,RCC_LR)

% [Slope2,inter2,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval2] = TractSF(FAR(1,:),MutotlistAP,MutotlistLR,UperctotAP,UperctotLR,FAtractlistAP,FAtractlistLR,AmplistAP,AmplistLR);

FAtractlistAP1 = cat(1,FAtractlistAP1,FAtractlistAP);
FAtractlistLR1 = cat(1,FAtractlistLR1,FAtractlistLR);
UperctotAP1 = cat(1,UperctotAP1,UperctotAP);
UperctotLR1 = cat(1,UperctotLR1,UperctotLR);
MutotlistAP1 = cat(1,MutotlistAP1,MutotlistAP);
MutotlistLR1 = cat(1,MutotlistLR1,MutotlistLR);
AmplistAP1 = cat(1,AmplistAP1,AmplistAP);
AmplistLR1 = cat(1,AmplistLR1,AmplistLR);

dir
cd /Volumes/CLJ-001-1/Group/drsmitty/ScanData/Multi/170222

cd Multiexcite_03; cd MRE; cd AP
load Quant_MUFS_AP.mat
cd ..; cd LR
load Quant_MUFS_LR.mat
cd ..
FAmin = .25;
FAmax = .85;
R_fa = zeros(13,43);
R_fa_AP = zeros(13,61);
R_fa_LR = zeros(13,61);

for Tract = 1:13
    Tract
    [Rfa,Rfa_AP,Rfa_LR,FAall] = FA_correlation(FAmin,FAmax,Tract,FAtractlistAP,FAtractlistLR,MutotlistAP,UperctotAP,MutotlistLR,UperctotLR,AmplistAP,AmplistLR);
    R_fa(Tract,1:length(Rfa)) = abs(Rfa);
    R_fa_AP(Tract,1:length(Rfa)) = abs(Rfa_AP);
    R_fa_LR(Tract,1:length(Rfa)) = abs(Rfa_LR);
    
    R(1,Tract) = max(abs(Rfa));
    R(2,Tract) = max(abs(Rfa_AP));
    R(3,Tract) = max(abs(Rfa_LR));    
    
    Rfaall = FAall(Rfa == max(abs(Rfa)));
    RfaLR = FAall(Rfa_LR == max(abs(Rfa_LR)));
    RfaAP = FAall(Rfa_AP == max(abs(Rfa_AP)));
    
    FAR(1,Tract) = Rfaall(1);
    FAR(2,Tract) = RfaAP(1);
    FAR(3,Tract) = RfaLR(1);    
end

RCC(1,:) = R_fa(1,:); RCC(2,:) = R_fa(6,:); RCC(3,:) = R_fa(13,:);
RCC_AP(1,:) = R_fa_AP(1,:); RCC(2,:) = R_fa_AP(6,:); RCC(3,:) = R_fa_AP(13,:);
RCC_LR(1,:) = R_fa_LR(1,:); RCC(2,:) = R_fa_LR(6,:); RCC(3,:) = R_fa_LR(13,:);
figure;plot(FAall,RCC)
figure;plot(FAall,RCC_AP)
figure;plot(FAall,RCC_LR)

% [Slope3,inter3,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval3] = TractSF(FAR(1,:),MutotlistAP,MutotlistLR,UperctotAP,UperctotLR,FAtractlistAP,FAtractlistLR,AmplistAP,AmplistLR);

FAtractlistAP1 = cat(1,FAtractlistAP1,FAtractlistAP);
FAtractlistLR1 = cat(1,FAtractlistLR1,FAtractlistLR);
UperctotAP1 = cat(1,UperctotAP1,UperctotAP);
UperctotLR1 = cat(1,UperctotLR1,UperctotLR);
MutotlistAP1 = cat(1,MutotlistAP1,MutotlistAP);
MutotlistLR1 = cat(1,MutotlistLR1,MutotlistLR);
AmplistAP1 = cat(1,AmplistAP1,AmplistAP);
AmplistLR1 = cat(1,AmplistLR1,AmplistLR);

cd /Volumes/CLJ-001-1/Group/drsmitty/ScanData/Multi/170816

cd Clj_Multiexcite_04_02; cd MRE; cd AP
load Quant_MUFS_AP.mat
cd ..; cd LR
load Quant_MUFS_LR.mat
cd ..
FAmin = .25;
FAmax = .85;
R_fa = zeros(13,43);
R_fa_AP = zeros(13,61);
R_fa_LR = zeros(13,61);

FAtractlistAP1 = cat(1,FAtractlistAP1,FAtractlistAP);
FAtractlistLR1 = cat(1,FAtractlistLR1,FAtractlistLR);
UperctotAP1 = cat(1,UperctotAP1,UperctotAP);
UperctotLR1 = cat(1,UperctotLR1,UperctotLR);
MutotlistAP1 = cat(1,MutotlistAP1,MutotlistAP);
MutotlistLR1 = cat(1,MutotlistLR1,MutotlistLR);
AmplistAP1 = cat(1,AmplistAP1,AmplistAP);
AmplistLR1 = cat(1,AmplistLR1,AmplistLR);

for Tract = 1:13
    Tract
    [Rfa,Rfa_AP,Rfa_LR,FAall] = FA_correlation(FAmin,FAmax,Tract,FAtractlistAP,FAtractlistLR,MutotlistAP,UperctotAP,MutotlistLR,UperctotLR,AmplistAP,AmplistLR);
    R_fa(Tract,1:length(Rfa)) = abs(Rfa);
    R_fa_AP(Tract,1:length(Rfa)) = abs(Rfa_AP);
    R_fa_LR(Tract,1:length(Rfa)) = abs(Rfa_LR);
    
    R(1,Tract) = max(abs(Rfa));
    R(2,Tract) = max(abs(Rfa_AP));
    R(3,Tract) = max(abs(Rfa_LR));    
    
    Rfaall = FAall(Rfa == max(abs(Rfa)));
    RfaLR = FAall(Rfa_LR == max(abs(Rfa_LR)));
    RfaAP = FAall(Rfa_AP == max(abs(Rfa_AP)));
    
    FAR(1,Tract) = Rfaall(1);
    FAR(2,Tract) = RfaAP(1);
    FAR(3,Tract) = RfaLR(1);    
end

RCC(1,:) = R_fa(1,:); RCC(2,:) = R_fa(6,:); RCC(3,:) = R_fa(13,:);
RCC_AP(1,:) = R_fa_AP(1,:); RCC(2,:) = R_fa_AP(6,:); RCC(3,:) = R_fa_AP(13,:);
RCC_LR(1,:) = R_fa_LR(1,:); RCC(2,:) = R_fa_LR(6,:); RCC(3,:) = R_fa_LR(13,:);
figure;plot(FAall,RCC)
figure;plot(FAall,RCC_AP)
figure;plot(FAall,RCC_LR)

[Slope,inter,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval] = TractSF(FAR(1,:),MutotlistAP1,MutotlistLR1,UperctotAP1,UperctotLR1,FAtractlistAP1,FAtractlistLR1,AmplistAP1,AmplistLR1);

% [Slope4,inter4,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval4] = TractSF(FAR(1,:),MutotlistAP,MutotlistLR,UperctotAP,UperctotLR,FAtractlistAP,FAtractlistLR,AmplistAP,AmplistLR);

MuS = inter(1:2,:);
MuF = MuS + 100*Slope(1:2,:);

% figure;bar(MuS)
% figure;bar(MuF)
% figure;bar(Slope(1:2,:)*100)
% figure;bar(Slope(3,:)*100)




