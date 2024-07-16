function [M_t,U_t_all,U_f_final,U_s_final,total_4D] = post_crss(U_s_all,U_f_all,Brfo_all,A,a)

Xmotion = Brfo_all(:,:,:,1);
Ymotion = Brfo_all(:,:,:,2);
Zmotion = Brfo_all(:,:,:,3);

N = Brfo_all;
S_x_s = length(Xmotion(:,1,1));
S_y_s = length(Ymotion(1,:,1));
S_z_s = length(Zmotion(1,1,:));

M_t = abs(U_f_all + U_s_all);

cd('/Users/drsmitty/Desktop/Epilepsy_ISMRM/ISMRM/Multiexcite_02/MRE/Tracts')

CC = load_nii('tractsCC.nii'); CC = CC.img; CC = double(CC>0.7);
CR = load_nii('tractsCR.nii'); CR = CR.img; CR = double(CR>0.7);
SLF = load_nii('tractsSLF.nii'); SLF = SLF.img; SLF = double(SLF>0.7);

total = double((SLF+CR+CC)>0.7);
total_4D(:,:,:,1) = total; total_4D(:,:,:,2) = total; total_4D(:,:,:,3) = total; 
tmp = load_nii('ref.nii');
tmp = tmp.img;
tmpflip1 = permute(tmp,[2,1,3]);
U_s_tot = zeros(size(U_s_all(:,:,:,1)));
U_f_tot = zeros(size(U_s_tot));
U_t_tot = zeros(size(U_s_tot));
U_t_all = total_4D.*M_t;
tic
for i = 1:S_x_s
    for j = 1:S_y_s
        for k = 1:S_z_s
            U_s_tot(i,j,k) = sqrt(U_s_all(i,j,k,1)^2+U_s_all(i,j,k,2)^2+U_s_all(i,j,k,3)^2);
            U_f_tot(i,j,k) = sqrt(U_f_all(i,j,k,1)^2+U_f_all(i,j,k,2)^2+U_f_all(i,j,k,3)^2);
            U_t_tot(i,j,k) = sqrt(U_t_all(i,j,k,1)^2+U_t_all(i,j,k,2)^2+U_t_all(i,j,k,3)^2);
        end
    end
end
toc
U_t_tot = U_t_tot.*total;

U_f_final = zeros(size(U_s_tot));
U_s_final = zeros(size(U_s_tot));
tic
for i = 1:S_x_s
    for j = 1:S_y_s
        for k = 1:S_z_s
            if U_t_tot(i,j,k) == 0
                U_s_final(i,j,k) = NaN;
                U_f_final(i,j,k) = NaN;
            else
                U_s_final(i,j,k) = U_s_tot(i,j,k)/U_t_tot(i,j,k);
                U_f_final(i,j,k) = U_f_tot(i,j,k)/U_t_tot(i,j,k);
            end
        end
    end
end
toc
total_4D(:,:,:,1) = total; total_4D(:,:,:,2) = total; total_4D(:,:,:,3) = total; 

sprintf('making figures')

figure;
imshow(abs(flip(permute(U_s_final(:,:,30),[2 1 3]),1))); caxis([0 1]); set(gcf,'colormap',hot);colorbar;
figure;
imshow(abs(flip(permute(U_f_final(:,:,30),[2 1 3]),1))); caxis([0 1]); set(gcf,'colormap',hot);colorbar;
figure;
imshow(abs(flip(permute(U_t_tot(:,:,30),[2 1 3]),1))); caxis([0 .75]); set(gcf,'colormap',hot);colorbar;
% figure;
% imshow(tmp(:,:,30)); caxis([0 8000]); set(gcf,'colormap',gray);colorbar;
% figure;
% imshow(total(:,:,30)); caxis([0 1]); set(gcf,'colormap',gray);colorbar;
% figure;
% imshow(CC(:,:,30)); caxis([0 1]); set(gcf,'colormap',gray);colorbar;
% figure;
% imshow(CR(:,:,30)); caxis([0 1]); set(gcf,'colormap',gray);colorbar;
% figure;
% imshow(SLF(:,:,30)); caxis([0 1]); set(gcf,'colormap',gray);colorbar;

% U_s_ROI = U_s_all.*total_4D;
% U_f_ROI = U_f_all.*total_4D;
% A = A.*total_4D;
% quiv(U_s_ROI,S_x_s,S_z_s,a,3)
% quiv(U_f_ROI,S_x_s,S_z_s,a,3)
% quiv(A,S_x_s,S_z_s,a,1)

function quiv(Output,xy,z,a,scaling)
tmp = load_nii('ref.nii');
tmp = tmp.img;
tmpflip1 = tmp;

nskip = a;
X = repmat((1:xy),[xy 1 z]);
Y = repmat((1:xy)',[1 xy z]);
Z = zeros(size(X));
for ii = 1:z
Z(:,:,ii) = ii;
end


UX = Output(:,:,:,1); UX = UX.*scaling; UX = UX(1:nskip:xy,1:nskip:xy,1:nskip:z); 
UY = Output(:,:,:,2); UY = UY.*scaling; UY = UY(1:nskip:xy,1:nskip:xy,1:nskip:z); 
UZ = Output(:,:,:,3); UZ = UZ.*scaling; UZ = UZ(1:nskip:xy,1:nskip:xy,1:nskip:z); 

S = size(UX);
S_x = S(1); S_y = S(2); S_z = S(3);
qx = zeros(S_x,S_y,S_z); qy = zeros(S_x,S_y,S_z); qz = zeros(S_x,S_y,S_z);

for i = 1:S_x
   for j = 1:S_y
        for k = 1:S_z          
            if abs(UX(i,j,k)) > abs(UY(i,j,k)) && abs(UX(i,j,k)) > abs(UZ(i,j,k))
               qx(i,j,k) = 1;               
            end  
            if abs(UY(i,j,k)) > abs(UX(i,j,k)) && abs(UY(i,j,k)) > abs(UZ(i,j,k))
               qy(i,j,k) = 1;         
            end 
            if abs(UZ(i,j,k)) > abs(UX(i,j,k)) && abs(UZ(i,j,k)) > abs(UY(i,j,k))
               qz(i,j,k) = 1;                         
            end 
        end
   end
end

% mask = ones(S_x,S_y,S_z); oneslice = ones(S_x,S_y);
% mask(:,:,1) = mask(:,:,1) - oneslice; mask(:,:,2) = mask(:,:,2) - oneslice;
% mask(:,:,S_z) = mask(:,:,S_z) - oneslice; mask(:,:,S_z-1) = mask(:,:,S_z-1) - oneslice;

XX = X(1:nskip:xy,1:nskip:xy,1:nskip:z);
YY = Y(1:nskip:xy,1:nskip:xy,1:nskip:z);
ZZ = Z(1:nskip:xy,1:nskip:xy,1:nskip:z);

UXX = UX.*qx; UXY = UY.*qx; UXZ = UZ.*qx;
UYX = UX.*qy; UYY = UY.*qy; UYZ = UZ.*qy;
UZX = UX.*qz; UZY = UY.*qz; UZZ = UZ.*qz;

% zeromask = ones(size(XX));
% for i = 1:S_x
%     for j = 1:S_y
%         for k = 1:S_z
%             if UXX(i,j,k) == 0 || UYX(i,j,k) == 0 || UZX(i,j,k) == 0 || UXY(i,j,k) == 0 || UXZ(i,j,k) == 0 || UYY(i,j,k) == 0 || UYZ(i,j,k) == 0 || UZZ(i,j,k) == 0 || UZY(i,j,k) == 0
%                 zeromask(i,j,k) = 0;
%             elseif isnan(UXX1(i,j,k)) == 0 || isnan(UYX1(i,j,k)) == 0 || isnan(UZX1(i,j,k)) == 0 || isnan(UXY1(i,j,k)) == 0 || isnan(UXZ1(i,j,k)) == 0 || isnan(UYY1(i,j,k)) == 0 || isnan(UYZ1(i,j,k)) == 0 || isnan(UZZ(i,j,k)) == 0 || isnan(UZY(i,j,k)) == 0
%                 zeromask(i,j,k) = 0;
%             end
%         end
%     end
% end
% 
% UZX(UZX==0) = nan;UXX(UXX==0) = nan;UYX(UYX==0) = nan;
% UZZ(UZZ==0) = nan;UXZ(UXZ==0) = nan;UYZ(UYZ==0) = nan;
% UZY(UZY==0) = nan;UXY(UXY==0) = nan;UYY(UYY==0) = nan;
% 
% figure;
% hold on
% quiver3(XX,YY,ZZ,UXX,UXY,UXZ,2,'r'); axis equal;
% quiver3(XX,YY,ZZ,UYX,UYY,UYZ,2,'color',[0 .65 0]); axis equal;
% quiver3(XX,YY,ZZ,UZX,UZY,UZZ,2,'b'); axis equal;
% hold off

figure;
hold on
quiver3(YY,XX,ZZ,UXX,UXY,UXZ,1,'r'); axis equal;
quiver3(YY,XX,ZZ,UYX,UYY,UYZ,1,'color',[0 .65 0]); axis equal;
quiver3(YY,XX,ZZ,UZX,UZY,UZZ,1,'b'); axis equal;
hold off


% figure;
% hold on
% im(tmpflip1(:,:,30));
% quiver(XX(:,:,30),YY(:,:,30),UXX(:,:,30),UXY(:,:,30),3,'r'); axis equal;
% quiver(XX(:,:,30),YY(:,:,30),UYX(:,:,30),UYY(:,:,30),3,'color',[0 .65 0]); axis equal;
% quiver(XX(:,:,30),YY(:,:,30),UZX(:,:,30),UZY(:,:,30),3,'b'); axis equal;
% hold off

figure;
hold on
im(tmpflip1(:,:,30));
quiver(YY(:,:,30),XX(:,:,30),UXX(:,:,30),UXY(:,:,30),1.5,'r'); axis equal;
quiver(YY(:,:,30),XX(:,:,30),UYX(:,:,30),UYY(:,:,30),1.5,'color',[0 .65 0]); axis equal;
quiver(YY(:,:,30),XX(:,:,30),UZX(:,:,30),UZY(:,:,30),1.5,'b'); axis equal;
hold off
