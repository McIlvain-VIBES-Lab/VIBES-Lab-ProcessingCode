Uperctot = cat(4,UperctotAP,UperctotLR);
Mutotlist = cat(3,MutotlistAP,MutotlistLR);
Thlist = cat(4,ThlistAP,ThlistLR);
vox = zeros(13,3,length(Mutotlist(1,1,:)));

for ti = 1:13;
    if ti == 3 || ti == 6 || ti == 13
        clear b_new Sv1 Fv1 AS AF MuS MuF v1 Mu MuT
        ti
        SF_thresh = 60;
        FA_thresh = 0.6;
        FAv = FAtractlistAP(:,ti)>FA_thresh;
        
        b = [];
        A = [];
        for dirs = 1:length(Mutotlist(1,1,:))
            
            clear Sv1 Fv1 AS AF MuS MuF v1 Mu MuT
            Sv1 = Uperctot(:,ti,1,dirs)>=SF_thresh; Fv1 = Uperctot(:,ti,2,dirs)>=SF_thresh;
            MuS = Mutotlist(logical(Sv1.*FAv),ti,dirs); MuF = Mutotlist(logical(Fv1.*FAv),ti,dirs);
            v1 = Sv1+Fv1;
            MuT = Mutotlist(:,ti,dirs); Mu = MuT(logical((MuT>0).*FAv.*v1));
            
            AS(1:length(MuS),1) = 1;
            AS(1:length(MuS),2) = cos(Thlist(logical(Sv1.*FAv),1,ti,dirs)).^2;
            AS(1:length(MuS),3) = 0;
            
            AF(1:length(MuF),1) = 1;
            AF(1:length(MuF),2) = cos(2*Thlist(logical(Fv1.*FAv),1,ti,dirs)).^2;
            AF(1:length(MuF),3) = sin(2*Thlist(logical(Fv1.*FAv),1,ti,dirs)).^2;
            
            Arep = cat(1,(AS.*repmat(squeeze(Uperctot(logical(Sv1.*FAv),ti,1,dirs)),[1 3])/100),(AF.*repmat(squeeze(Uperctot(logical(Fv1.*FAv),ti,2,dirs)),[1 3]))/100);
            A = cat(1,A,Arep);
            
            b = cat(1,b,Mu);
            vox(ti,dirs) = length(Mu);
            
        end
        size(b)
        size(A)
        
        x = A\b;
        b_new = A*x;
        
        x_soln(:,ti) = x;
        r2(ti,1) = rsquare(b,b_new);
        clear A b x
    end
end

