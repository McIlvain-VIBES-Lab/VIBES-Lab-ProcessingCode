function [Rfa,Rfa_AP,Rfa_LR,FAall] = FA_correlation(FAmin,FAmax,T,FAtractlistAP,FAtractlistLR,MutotlistAP,UperctotAP,MutotlistLR,UperctotLR,AmplistAP,AmplistLR)

ii = 0;
Rfa_AP = [];
Rfa_LR = [];
Rfa = [];

for fa = FAmin:.01:FAmax
    ii = ii+1;
    FA_thresh = fa;
    FAvAP = FAtractlistAP(:,T)>FA_thresh;
    FAvLR = FAtractlistLR(:,T)>FA_thresh;
    tmp = MutotlistAP(logical(FAvAP),T);
    tmp2 = MutotlistLR(logical(FAvAP),T);
    if sum(tmp>0) >= 25
        if sum(tmp2>0) >= 25
            MutractAP = MutotlistAP(logical(FAvAP),T);
            Sperctot1AP = UperctotAP(logical(FAvAP),T,1);
            Sperctot2AP = UperctotAP(logical(FAvAP),T,3);
            Fperctot1AP = UperctotAP(logical(FAvAP),T,2);
            Fperctot2AP = UperctotAP(logical(FAvAP),T,4);
            Amp1AP = AmplistAP(logical(FAvAP),T,1);
            Amp2AP = AmplistAP(logical(FAvAP),T,2);
            
            S1 = Amp1AP.*Sperctot1AP; S2 = Amp2AP.*Sperctot2AP;
            F1 = Amp1AP.*Fperctot1AP; F2 = Amp2AP.*Fperctot2AP;
            APS = S1+S2; APF = F1+F2;
            
            MutractLR = MutotlistLR(logical(FAvLR),T);
            Sperctot1LR = UperctotLR(logical(FAvLR),T,1);
            Sperctot2LR = UperctotLR(logical(FAvLR),T,3);
            Fperctot1LR = UperctotLR(logical(FAvLR),T,2);
            Fperctot2LR = UperctotLR(logical(FAvLR),T,4);
            Amp1LR = AmplistLR(logical(FAvLR),T,1);
            Amp2LR = AmplistLR(logical(FAvLR),T,2);
            
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
            
            [RAP] = corr(APF,MutractAP);
            [RLR] = corr(LRF,MutractLR);
            [R] = corr(F,Mutract);
            Rfa_AP(ii) = abs(RAP);
            Rfa_LR(ii) = abs(RLR);
            Rfa(ii) = abs(R);
        else
            
        end
    else
        Rfa_AP(ii) = 0;
        Rfa_LR(ii) = 0;
        Rfa(ii) = 0;
    end
    FAall(ii) = fa;
end



% figure;plot(FAall,Rfa)
% figure;plot(FAall,Rfa_AP)
% figure;plot(FAall,Rfa_LR)
