% two excitation, slow/fast discrimination on V1
for ti = 1:13;
    if ti ~= 4 && ti ~= 5
        clear MuAP MuLR MuTAP MuTLR b ASAP ASLR AFAP AFLR AAP ALR 
        ti
  
        MuTAP = MutotlistAP(:,ti);
        MuAP = MuTAP(MuTAP>0);
        
        MuTLR = MutotlistLR(:,ti);
        MuLR = MuTLR(MuTLR>0);
        
        %MuTSI= MutotlistSI(:,ti);
        %MuSI = MuTSI(MuTSI>0);
        
        b = [MuAP;MuLR];
        
        ASAP(1:length(MuAP),1) = 1;
        ASAP(1:length(MuAP),2) = cos(ThlistAP(1:length(MuAP),1,ti)).^2;
        ASAP(1:length(MuAP),3) = 0;
        
        AFAP(1:length(MuLR),1) = 1;
        AFAP(1:length(MuAP),2) = cos(2*ThlistAP(1:length(MuAP),1,ti)).^2;
        AFAP(1:length(MuAP),3) = sin(2*ThlistAP(1:length(MuAP),1,ti)).^2;
        
        ASLR(1:length(MuLR),1) = 1;
        ASLR(1:length(MuLR),2) = cos(ThlistLR(1:length(MuLR),1,ti)).^2;
        ASLR(1:length(MuLR),3) = 0;
        
        AFLR(1:length(MuLR),1) = 1;
        AFLR(1:length(MuLR),2) = cos(2*ThlistLR(1:length(MuLR),1,ti)).^2;
        AFLR(1:length(MuLR),3) = sin(2*ThlistLR(1:length(MuLR),1,ti)).^2;
        
        %ASSI(1:length(MuSI),1) = 1;
        %ASSI(1:length(MuSI),2) = cos(ThlistSI(1:length(MuSI),1,ti)).^2;
        %ASSI(1:length(MuSI),3) = 0;
        
        %AFSI(1:length(MuSI),1) = 1;
        %AFSI(1:length(MuSI),2) = cos(2*ThlistLR(1:length(MuSI),1,ti)).^2;
        %AFSI(1:length(MuSI),3) = sin(2*ThlistLR(1:length(MuSI),1,ti)).^2;
        
        AAP = (ASAP.*repmat(UperctotAP(1:length(MuAP),ti,1),[1 3])/100)+(AFAP.*repmat(UperctotAP(1:length(MuAP),ti,2),[1 3])/100);
        ALR = (ASLR.*repmat(UperctotLR(1:length(MuLR),ti,1),[1 3])/100)+(AFLR.*repmat(UperctotLR(1:length(MuLR),ti,2),[1 3])/100);
        %ASI = (ASSI.*repmat(UperctotSI(1:length(MuSI),ti,1),[1 3])/100)+(AFSI.*repmat(UperctotSI(1:length(MuSI),ti,2),[1 3])/100);
        
        A = cat(1,AAP,ALR);
        size(b)
        size(A)
        
        x = A\b;
        b_new = A*x;
        
        x_soln(:,ti) = x;
        r2(ti,1) = rsquare(b,b_new);
        clear A b b_new x MuY MuT
    end
end

