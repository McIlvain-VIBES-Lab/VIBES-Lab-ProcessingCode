function [Slope,inter,APMu,LRMu,Mu,APslow,APfast,LRslow,LRfast,slow,fast,Rval] = TractSF(FA_thresh,MutotlistAP,MutotlistLR,UperctotAP,UperctotLR,FAtractlistAP,FAtractlistLR,AmplistAP,AmplistLR)
% cd AP
% load Quant_MUFS_AP.mat
% cd ..; cd LR
% load Quant_MUFS_LR.mat
% cd ..

Toptions = [1 6 13];
APMu = zeros(length(MutotlistAP),length(Toptions));
LRMu = zeros(length(MutotlistAP),length(Toptions));
Mu = zeros(length(MutotlistAP),length(Toptions));
APslow = zeros(length(MutotlistAP),length(Toptions));
APfast = zeros(length(MutotlistAP),length(Toptions));
LRslow = zeros(length(MutotlistAP),length(Toptions));
LRfast = zeros(length(MutotlistAP),length(Toptions));
slow = zeros(length(MutotlistAP),length(Toptions));
fast = zeros(length(MutotlistAP),length(Toptions));
Pval = []; Rval = [];

for ii = 1:length(Toptions)
    
    T = Toptions(ii);   
    
    %% Combined
    FAvAPc = FAtractlistAP(:,T)>.6;
    FAvLRc = FAtractlistLR(:,T)>.6;
%     FAvAPc = FAtractlistAP(T,:)>.6;
%     FAvLRc = FAtractlistLR(T,:)>.6;

    
    MutractAP = MutotlistAP(logical(FAvAPc),T);
    Sperctot1AP = UperctotAP(logical(FAvAPc),T,1);
    Sperctot2AP = UperctotAP(logical(FAvAPc),T,3);
    Fperctot1AP = UperctotAP(logical(FAvAPc),T,2);
    Fperctot2AP = UperctotAP(logical(FAvAPc),T,4);
    Amp1AP = AmplistAP(logical(FAvAPc),T,1);
    Amp2AP = AmplistAP(logical(FAvAPc),T,2);
    
%     MutractAP = MutotlistAP(logical(FAvAPc),T);
%     Sperctot1AP = UperctotAP(logical(FAvAPc),1);
%     Sperctot2AP = UperctotAP(logical(FAvAPc),3);
%     Fperctot1AP = UperctotAP(logical(FAvAPc),2);
%     Fperctot2AP = UperctotAP(logical(FAvAPc),4);
%     Amp1AP = AmplistAP(logical(FAvAPc),1);
%     Amp2AP = AmplistAP(logical(FAvAPc),2);
%     
    S1 = Amp1AP.*Sperctot1AP; S2 = Amp2AP.*Sperctot2AP;
    F1 = Amp1AP.*Fperctot1AP; F2 = Amp2AP.*Fperctot2AP;
    APS = S1+S2; APF = F1+F2;
    
    MutractLR = MutotlistLR(logical(FAvLRc),T);
    Sperctot1LR = UperctotLR(logical(FAvLRc),T,1);
    Sperctot2LR = UperctotLR(logical(FAvLRc),T,3);
    Fperctot1LR = UperctotLR(logical(FAvLRc),T,2);
    Fperctot2LR = UperctotLR(logical(FAvLRc),T,4);
    Amp1LR = AmplistLR(logical(FAvLRc),T,1);
    Amp2LR = AmplistLR(logical(FAvLRc),T,2);

%     MutractLR = MutotlistLR(logical(FAvLRc),T);
%     Sperctot1LR = UperctotLR(logical(FAvLRc),1);
%     Sperctot2LR = UperctotLR(logical(FAvLRc),3);
%     Fperctot1LR = UperctotLR(logical(FAvLRc),2);
%     Fperctot2LR = UperctotLR(logical(FAvLRc),4);
%     Amp1LR = AmplistLR(logical(FAvLRc),1);
%     Amp2LR = AmplistLR(logical(FAvLRc),2);
    
    S1 = Amp1LR.*Sperctot1LR; S2 = Amp2LR.*Sperctot2LR;
    F1 = Amp1LR.*Fperctot1LR; F2 = Amp2LR.*Fperctot2LR;
    LRS = S1+S2; LRF = F1+F2;
    
    Amp1 = cat(1,Amp1AP,Amp1LR);
    Amp2 = cat(1,Amp2AP,Amp2LR);
    Mutract = cat(1,MutractAP,MutractLR);
    Sperctot1 = cat(1,Sperctot1AP,Sperctot1LR);
    Sperctot2 = cat(1,Sperctot2AP,Sperctot2LR);
    Fperctot1 = cat(1,Fperctot1AP,Fperctot1LR);
    Fperctot2 = cat(1,Fperctot2AP,Fperctot2LR);
    
    S1 = Amp1.*Sperctot1; S2 = Amp2.*Sperctot2;
    F1 = Amp1.*Fperctot1; F2 = Amp2.*Fperctot2;
    S = S1+S2; F = F1+F2;
    
    %% Reassign FA tracts
    FAvAP = FAtractlistAP(:,T)>.6;
    FAvLR = FAtractlistLR(:,T)>.6;
%     FAvAP = FAtractlistAP(T,:)>.6;
%     FAvLR = FAtractlistLR(T,:)>.6;
    
    %% AP
    MutractAP = MutotlistAP(logical(FAvAPc),T);
    Sperctot1AP = UperctotAP(logical(FAvAPc),T,1);
    Sperctot2AP = UperctotAP(logical(FAvAPc),T,3);
    Fperctot1AP = UperctotAP(logical(FAvAPc),T,2);
    Fperctot2AP = UperctotAP(logical(FAvAPc),T,4);
    Amp1AP = AmplistAP(logical(FAvAPc),T,1);
    Amp2AP = AmplistAP(logical(FAvAPc),T,2);
%     
%     MutractAP = MutotlistAP(logical(FAvAPc),T);
%     Sperctot1AP = UperctotAP(logical(FAvAPc),1);
%     Sperctot2AP = UperctotAP(logical(FAvAPc),3);
%     Fperctot1AP = UperctotAP(logical(FAvAPc),2);
%     Fperctot2AP = UperctotAP(logical(FAvAPc),4);
%     Amp1AP = AmplistAP(logical(FAvAPc),1);
%     Amp2AP = AmplistAP(logical(FAvAPc),2);
    
    
    S1 = Amp1AP.*Sperctot1AP; S2 = Amp2AP.*Sperctot2AP;
    F1 = Amp1AP.*Fperctot1AP; F2 = Amp2AP.*Fperctot2AP;
    APS = S1+S2; APF = F1+F2;
    
    %% LR
    MutractLR = MutotlistLR(logical(FAvLRc),T);
    Sperctot1LR = UperctotLR(logical(FAvLRc),T,1);
    Sperctot2LR = UperctotLR(logical(FAvLRc),T,3);
    Fperctot1LR = UperctotLR(logical(FAvLRc),T,2);
    Fperctot2LR = UperctotLR(logical(FAvLRc),T,4);
    Amp1LR = AmplistLR(logical(FAvLRc),T,1);
    Amp2LR = AmplistLR(logical(FAvLRc),T,2);
% 
%     MutractLR = MutotlistLR(logical(FAvLRc),T);
%     Sperctot1LR = UperctotLR(logical(FAvLRc),1);
%     Sperctot2LR = UperctotLR(logical(FAvLRc),3);
%     Fperctot1LR = UperctotLR(logical(FAvLRc),2);
%     Fperctot2LR = UperctotLR(logical(FAvLRc),4);
%     Amp1LR = AmplistLR(logical(FAvLRc),1);
%     Amp2LR = AmplistLR(logical(FAvLRc),2);
    
    
    S1 = Amp1LR.*Sperctot1LR; S2 = Amp2LR.*Sperctot2LR;
    F1 = Amp1LR.*Fperctot1LR; F2 = Amp2LR.*Fperctot2LR;
    LRS = S1+S2; LRF = F1+F2;
       
    %% Figures
%     figure; axis([20 80 2000 4500])
%     hold on
%     scatter(APS,MutractAP,'red'); 
%     scatter(LRS,MutractLR,'blue');
%     hold off
    
    figure; axis([20 80 2000 4500])
    hold on
    scatter(APF,MutractAP,'red'); 
    scatter(LRF,MutractLR,'blue');
    hold off
 
    %% Correlation Analysis
    [RAP,PAP] = corr(APF,MutractAP);
    Pfit = polyfit(APF,MutractAP,1);
    PolyAP = Pfit(1);
    crossAP = Pfit(2);
    %R2 = RAP^2;
    %figure;scatter(LRS,MutractLR); figure;scatter(LRF,MutractLR);
    [RLR,PLR] = corr(LRF,MutractLR);
    Pfit = polyfit(LRF,MutractLR,1);
    PolyLR = Pfit(1);
    crossLR = Pfit(2);
    %R2 = RLR^2;
    %figure;scatter(S,Mutract); figure;scatter(F,Mutract); 
    [R,P] = corr(F,Mutract);
    Pfit = polyfit(F,Mutract,1);
    Poly = Pfit(1);
    cross = Pfit(2);
    %R2 = R^2;
    
    APMu(1:length(MutractAP),ii) = MutractAP; 
    LRMu(1:length(MutractLR),ii) = MutractLR;
    Mu(1:length(Mutract),ii) = Mutract;
    APslow(1:length(APS),ii) = APS; APfast(1:length(APF),ii) = APF;
    LRslow(1:length(LRS),ii) = LRS; LRfast(1:length(LRF),ii) = LRF;
    slow(1:length(S),ii) = S; fast(1:length(F),ii) = F;
    
    Rval(1,ii) = RAP; Rval(2,ii) = RLR; Rval(3,ii) = R;
    Pval(1,ii) = PAP; Pval(2,ii) = PLR; Pval(3,ii) = P;
    
    Slope(1,ii) = PolyAP;
    Slope(2,ii) = PolyLR;
    Slope(3,ii) = Poly;
    
    inter(1,ii) = crossAP;
    inter(2,ii) = crossLR;
    inter(3,ii) = cross;
    
end