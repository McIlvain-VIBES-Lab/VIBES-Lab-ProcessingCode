function [x_soln, crit_out, x_soln1] = combo_2waves_LDI_crit(directions,threshold,waves)
%% Weighting the A with the amplitudes of the primary and secondary waves
% Includes Slow/Fast Thresholding & FA Thresholding option
%       threshold = 1 for on = 0 for off
% Includes option for # of directions
%       directions = 1, 2, or 3
% Written By DRS
% Finished 1/10/2018

% cd Skeleton
load FA.mat
% cd ..;
cd MRE; cd AP
load Quant_MUFS_APcrit.mat
cd ..

if directions == 1
    Mutotlist = MutotlistAP;
    Thlist = ThlistAP;
    crit = critAP;
elseif directions == 2
    cd LR; load Quant_MUFS_LRcrit.mat; cd ..
    Mutotlist = cat(4,MutotlistAP,MutotlistLR);
    Thlist = cat(4,ThlistAP,ThlistLR);
    crit = cat(2,critAP,critLR);
end
cd ..

vox = zeros(13,3,length(Mutotlist(1,1,:)));
ball = [];
% 
% cd Skeleton
% load FA.mat
% cd ..; cd MRE; cd AP
% load Quant_MUFS_APnew414.mat
% cd ..
% 
% if directions == 1
%     Uperctot = UperctotAP; Mutotlist = MutotlistAP;
%     Amptotlist = AmplistAP; Thlist = ThlistAP;
% elseif directions == 2
%     cd LR; load Quant_MUFS_LRnew414.mat; cd ..
%     Uperctot = cat(4,UperctotAP,UperctotLR); Mutotlist = cat(4,MutotlistAP,MutotlistLR);
%     Amptotlist = cat(4,AmplistAP,AmplistLR); Thlist = cat(4,ThlistAP,ThlistLR);
% end
% cd ..

% vox = zeros(13,3,length(Mutotlist(1,1,:)));
% ball = [];

for ti = 1:4;
    Adirs = [];
    b = [];
    if threshold == 1
        FA_thresh = .5;
    else
        FA_thresh = 0;
    end
    FAv = FAtractlistAP(:,ti)>FA_thresh;
    critall = sum(crit(:,:,ti),2)>0;
    tmp = critall.*FAv;
    % tmp = sum(critall);
    critall = repmat(tmp,[1 2*directions*waves]);
    
    clear Atot

        for wv = 1:waves
            for dirs = 1:length(Mutotlist(1,1,1,:))
                clear AS AF Mus Muf
                
                Musx = Mutotlist(:,ti,2*wv-1,dirs); Mufx = Mutotlist(:,ti,2*wv-1,dirs);
                Mus = Musx(logical((Musx>0).*FAv.*tmp)); Muf = Mufx(logical((Mufx>0).*FAv.*tmp));
                
                AS(1:length(Mus),1,1) = 1;
                AS(1:length(Mus),1,2) = cos(Thlist(logical(FAv.*tmp),wv,ti,dirs)).^2;
                AS(1:length(Mus),1,3) = 0;
                
                AF(1:length(Muf),1,1) = 1;
                AF(1:length(Muf),1,2) = cos(2*Thlist(logical(FAv.*tmp),wv,ti,dirs)).^2;
                AF(1:length(Muf),1,3) = sin(2*Thlist(logical(FAv.*tmp),wv,ti,dirs)).^2;
                
                Arep = cat(2,AS,AF);
                Adirs = cat(2,Adirs,Arep);
                
                Mu = cat(2,Mus,Muf);
                b = cat(2,b,Mu);
                vox(ti,dirs) = length(Mu);
            end
        end
        A = Adirs;
%         x = A\b
%         b_new = A*x;
        
 %%       this is the way to do each one individually
%         x = zeros(3,1);
%         if directions == 2
%             np = 8;
%         else
%             np = 4;
%         end
%         for tt = 1:(length(A)/np)
%             x = A(tt:(length(A)/np):end,:)\b(tt:(length(A)/np):end,1);
%             b_new = A(tt:(length(A)/np):end,:)*x;
%             x_soln(tt,:,ti) = x';
%             r2_sol(tt,1,ti) = rsquare(b(tt:(length(A)/np):end,1),b_new);
%         end
     
%         r2(1,ti) = rsquare(b,b_new);
%         ball(1:length(b),ti,1) = b;
%         ball(1:length(b_new),ti,2) = b_new;

%         critall = sum(crit(:,:,ti),2)>0;
%         tmp = critall.*FAv;
%         % tmp = sum(critall);
%         critall = repmat(tmp,[1 2*directions*waves]);
        
        criteria = [];
        for ii = 1:length(critall)
            if tmp(ii) == 1
                criteria = cat(1,criteria,crit(ii,1:2*waves*directions,ti));
            end     
        end  
        critrep = repmat(criteria,[1 1 3]);
        
        for tt = 1:size(A,1)
            b_new = squeeze(b(tt,:))'; % b(tt,logical(crits(tt,:)))
            A_new = squeeze(A(tt,:,:).*critrep(tt,:,:));
            x_soln(tt,:,ti) = A_new\b_new;
            x_soln1(tt,:,ti) = lsqnonneg(A_new,b_new);
            crit_out(tt,1,ti) = sum(critrep(tt,:,1),2);
        end
        
    clear A x
    %end
end

% x_soln
% r2
