for ti = 1:13;
    if ti ~= 4 && ti ~= 5
        clear b x A b_new MuLR MuAP MuSI Sv1_LR Fv1_LR Sv1_AP Fv1_AP AS_LR AF_LR AS_AP AF_AP MuS_LR MuF_LR MuS_AP MuF_AP Sv1_SI Fv1_SI AS_SI AF_SI MuS_SI MuF_SI v1_AP v1_LR v1_SI
        ti
        SF_thresh = 60;
        FA_thresh = 0.7;
        
        FAv = FAtractlistAP(:,ti)>FA_thresh;
        
        Sv1_AP = UperctotAP(:,ti,1)>=SF_thresh;
        Fv1_AP = UperctotAP(:,ti,2)>=SF_thresh;
        MuS_AP = MutotlistAP(logical(Sv1_AP.*FAv),ti);
        MuF_AP = MutotlistAP(logical(Fv1_AP.*FAv),ti);
        v1_AP = Sv1_AP+Fv1_AP;
        MuTAP = MutotlistAP(:,ti);
        MuAP = MuTAP(logical((MuTAP>0).*FAv.*v1_AP));
        
        Sv1_LR = UperctotLR(:,ti,1)>=SF_thresh;
        Fv1_LR = UperctotLR(:,ti,2)>=SF_thresh;
        MuS_LR = MutotlistLR(logical(Sv1_LR.*FAv),ti);
        MuF_LR = MutotlistLR(logical(Fv1_LR.*FAv),ti);
        v1_LR = Sv1_LR+Fv1_LR;
        MuTLR = MutotlistLR(:,ti);
        MuLR = MuTLR(logical((MuTLR>0).*FAv.*v1_LR));
        
%         Sv1_SI = UperctotSI(:,ti,1)>=SF_thresh;
%         Fv1_SI = UperctotSI(:,ti,2)>=SF_thresh;
%         MuS_SI = MutotlistSI(logical(Sv1_SI.*FAv),ti);
%         MuF_SI = MutotlistSI(logical(Fv1_SI.*FAv),ti);
%         v1_SI = Sv1_SI+Fv1_SI;
%         MuTSI= MutotlistSI(:,ti);
%         MuSI = MuTSI(logical((MuTSI>0).*FAv.*v1_SI));
        
%         b = cat(1,MuAP,MuLR,MuSI);
        b = cat(1,MuAP,MuLR);
        for dirs = 1:2
            Munam = Munames(dirs,:);
            voxtot = length(Munam);
            vox_num = [];
        end
        size(b)
        
        AS_AP(1:length(MuS_AP),1) = 1;
        AS_AP(1:length(MuS_AP),2) = cos(ThlistAP(logical(Sv1_AP.*FAv),1,ti)).^2;
        AS_AP(1:length(MuS_AP),3) = 0;
        
        AF_AP(1:length(MuF_AP),1) = 1;
        AF_AP(1:length(MuF_AP),2) = cos(2*ThlistAP(logical(Fv1_AP.*FAv),1,ti)).^2;
        AF_AP(1:length(MuF_AP),3) = sin(2*ThlistAP(logical(Fv1_AP.*FAv),1,ti)).^2;
        
        AS_LR(1:length(MuS_LR),1) = 1;
        AS_LR(1:length(MuS_LR),2) = cos(ThlistLR(logical(Sv1_LR.*FAv),1,ti)).^2;
        AS_LR(1:length(MuS_LR),3) = 0;
        
        AF_LR(1:length(MuF_LR),1) = 1;
        AF_LR(1:length(MuF_LR),2) = cos(2*ThlistLR(logical(Fv1_LR.*FAv),1,ti)).^2;
        AF_LR(1:length(MuF_LR),3) = sin(2*ThlistLR(logical(Fv1_LR.*FAv),1,ti)).^2;
        
%         AS_SI(1:length(MuS_SI),1) = 1;
%         AS_SI(1:length(MuS_SI),2) = cos(ThlistSI(logical(Sv1_SI.*FAv),1,ti)).^2;
%         AS_SI(1:length(MuS_SI),3) = 0;
%         
%         AF_SI(1:length(MuF_SI),1) = 1;
%         AF_SI(1:length(MuF_SI),2) = cos(2*ThlistSI(logical(Fv1_SI.*FAv),1,ti)).^2;
%         AF_SI(1:length(MuF_SI),3) = sin(2*ThlistSI(logical(Fv1_SI.*FAv),1,ti)).^2;
        
        AAP = cat(1,(AS_AP.*repmat(UperctotAP(logical(Sv1_AP.*FAv),ti,1),[1 3])/100),(AF_AP.*repmat(UperctotAP(logical(Fv1_AP.*FAv),ti,2),[1 3])/100));
        ALR = cat(1,(AS_LR.*repmat(UperctotLR(logical(Sv1_LR.*FAv),ti,1),[1 3])/100),(AF_LR.*repmat(UperctotLR(logical(Fv1_LR.*FAv),ti,2),[1 3])/100));
%         ASI = cat(1,(AS_SI.*repmat(UperctotSI(logical(Sv1_SI.*FAv),ti,1),[1 3])/100),(AF_SI.*repmat(UperctotSI(logical(Fv1_SI.*FAv),ti,2),[1 3])/100));
        
%         A = cat(1,AAP,ALR,ASI);
        A = cat(1,AAP,ALR);
        size(A)
        
        x = A\b;
        b_new = A*x;
        
        x_soln(:,ti) = x;
        r2(ti,1) = rsquare(b,b_new);
    end
end