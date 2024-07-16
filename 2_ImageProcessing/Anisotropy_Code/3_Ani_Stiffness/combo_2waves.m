function combo_2waves(directions,threshold,waves)
%% Weighting the A with the amplitudes of the primary and secondary waves
% Includes Slow/Fast Thresholding & FA Thresholding option
%       threshold = 1 for on = 0 for off
% Includes option for # of directions
%       directions = 1, 2, or 3
% Written By DRS
% Finished 1/10/2018

load FA.mat
cd MRE; cd AP
load Quant_MUFS_AP.mat
cd ..

if directions == 1
    Uperctot = UperctotAP; Mutotlist = MutotlistAP;
    Amptotlist = AmplistAP; Thlist = ThlistAP;
elseif directions == 2
    cd LR; load Quant_MUFS_LR.mat; cd ..
    Uperctot = cat(4,UperctotAP,UperctotLR); Mutotlist = cat(3,MutotlistAP,MutotlistLR);
    Amptotlist = cat(4,AmplistAP,AmplistLR); Thlist = cat(4,ThlistAP,ThlistLR);
elseif directions == 3
    cd LR; load Quant_MUFS_LR.mat; cd ..
    cd SI; load Quant_MUFS_SI.mat; cd ..
    Uperctot = cat(4,UperctotAP,UperctotLR,UperctotSI); Mutotlist = cat(3,MutotlistAP,MutotlistLR,MutotlistSI);
    Amptotlist = cat(4,AmplistAP,AmplistLR,AmplistSI); Thlist = cat(4,ThlistAP,ThlistLR,ThlistSI);
end
cd ..

vox = zeros(13,3,length(Mutotlist(1,1,:)));
ball = [];

for ti = 1:13;
    %if ti == 3 || ti == 6 || ti == 13
    %if ti >= 7
    Amp = [];
    Adirs = [];
    b = [];
    clear Atot
        if threshold == 1
            SF_thresh = 50;
            FA_thresh = .5;
        else 
            SF_thresh = 50;
            FA_thresh = 0;
        end
        FAv = FAtractlistAP(:,ti)>FA_thresh;
        for wv = 1:waves
            for dirs = 1:length(Mutotlist(1,1,:))
                clear Sv Fv AS AF MuS MuF v Mu MuT
                
                Sv = Uperctot(:,ti,(2*wv-1),dirs)>=SF_thresh; Fv = Uperctot(:,ti,(2*wv),dirs)>=SF_thresh;
                MuS = Mutotlist(logical(Sv.*FAv),ti,dirs); MuF = Mutotlist(logical(Fv.*FAv),ti,dirs);
                v = (Sv+Fv)>0;
                MuT = Mutotlist(:,ti,dirs); Mu = MuT(logical((MuT>0).*FAv.*v));
                
                AS(1:length(MuS),1) = 1;
                AS(1:length(MuS),2) = cos(Thlist(logical(Sv.*FAv),wv,ti,dirs)).^2;
                AS(1:length(MuS),3) = 0;
                
                AF(1:length(MuF),1) = 1;
                AF(1:length(MuF),2) = cos(2*Thlist(logical(Fv.*FAv),wv,ti,dirs)).^2;
                AF(1:length(MuF),3) = sin(2*Thlist(logical(Fv.*FAv),wv,ti,dirs)).^2;
                
                AmpS = Amptotlist(logical(Sv.*FAv),ti,wv,dirs);
                AmpF = Amptotlist(logical(Fv.*FAv),ti,wv,dirs);
                
                Arep = cat(1,(AS.*repmat(squeeze(Uperctot(logical(Sv.*FAv),ti,(2*wv-1),dirs)),[1 3])/100),(AF.*repmat(squeeze(Uperctot(logical(Fv.*FAv),ti,(2*wv),dirs)),[1 3]))/100);
                Adirs = cat(1,Adirs,Arep);
                
                Amptot = cat(1,AmpS,AmpF); 
                Amp = cat(1,Amp,Amptot);
                
                b = cat(1,b,Mu);
                vox(ti,dirs) = length(Mu);
            end
            if waves == 1
                Atot(1:length(Adirs),:,wv) = Adirs;
            else 
                Atot(1:length(Adirs),:,wv) = Adirs.*(repmat(Amp,[1 3]));
            end
        end
        
        %size(Atot)
        A = sum(Atot,3);
        %size(b)
        %size(A)
        % x = A\b;
%         x = zeros(size(A));
%         for tt = 1:length(Atot)
%             x(tt,:) = A(tt,:)\b(tt);
%         end
%         b_new = A*x;
        
        % this is the way to do each one individually
        x = zeros(3,1);
        for tt = 1:length(Atot)
            x = A(tt,:)\b(tt,1);
            b_new(tt,1) = A(tt,:)*x;
            x_soln(tt,:) = x';
        end
  
        
        
        r2(1,ti) = rsquare(b,b_new);
        ball(1:length(b),ti,1) = b;
        ball(1:length(b_new),ti,2) = b_new;
        clear A x
    %end
end

x_soln
r2
