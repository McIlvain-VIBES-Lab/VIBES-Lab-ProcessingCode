function [m_s,m_f] = crss(N,A,a)

xy = size(N,1); z = size(N,3);
m_f = zeros(xy,xy,z,3); m_s = zeros(xy,xy,z,3);

for i = 1:xy
    for j = 1:xy
        for k = 1:z
           tmp = N(i,j,k,:);
           tmp = squeeze(tmp);
           tmp1 = A(i,j,k,:);
           tmp1 = squeeze(tmp1);
           m_s(i,j,k,:) = cross(tmp,tmp1);
           tmp2 = m_s(i,j,k,:);
           tmp2 = squeeze(tmp2);
           m_f(i,j,k,:) = cross(tmp,tmp2);
        end
    end
end

quiv(m_s,xy,z,a,1)
quiv(m_f,xy,z,a,1)

function quiv(Output,xy,z,a,scaling)

cd('/Volumes/CLJ-001-2/Group/drsmitty/ScanData/Multi/170222/Multiexcite_02/MRE/Tracts')

CC = load_nii('tractsCC.nii'); CC = CC.img;
CR = load_nii('tractsCR.nii'); CR = CR.img;
SLF = load_nii('tractsSLF.nii'); SLF = SLF.img;
total = SLF+CR+CC;
tmp = load_nii('ref.nii');
tmp = tmp.img;
tmpflip1 = permute(tmp,[2,1,3]);

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

onem = ones(S_x,S_y,S_z); zerom = zeros(S_x,S_y,S_z);

%for i = 3:S_x-2
%    for j = 3:S_y-2
%        for k = 3:S_z-2
%            if UX(i+1,j,k) == 0 || UX(i-1,j,k) == 0 || UY(i+1,j,k) == 0 || UY(i-1,j,k) == 0 || UZ(i+1,j,k) == 0 || UZ(i-1,j,k) == 0 
%                zerom(i,j,k) = 1;
%            elseif UX(i,j+1,k) == 0 || UX(i,j-1,k) == 0 || UY(i,j+1,k) == 0 || UY(i,j-1,k) == 0 || UZ(i,j+1,k) == 0 || UZ(i,j-1,k) == 0 
%                zerom(i,j,k) = 1;
%            elseif UX(i,j,k+1) == 0 || UX(i,j,k-1) == 0 || UY(i,j,k+1) == 0 || UY(i,j,k-1) == 0 || UZ(i,j,k+1) == 0 || UZ(i,j,k-1) == 0 
%                zerom(i,j,k) = 1;
%            end
%        end
%    end
%end

oneslice = ones(S_x,S_y); mask = onem; %- zerom;
mask(:,:,1) = mask(:,:,1) - oneslice; mask(:,:,2) = mask(:,:,2) - oneslice;
mask(:,:,S_z) = mask(:,:,S_z) - oneslice; mask(:,:,S_z-1) = mask(:,:,S_z-1) - oneslice;

XX = X(1:nskip:xy,1:nskip:xy,1:nskip:z);
YY = Y(1:nskip:xy,1:nskip:xy,1:nskip:z);
ZZ = Z(1:nskip:xy,1:nskip:xy,1:nskip:z);
ROIs = [total; CC; CR; SLF];

for r = 1:4
    ROI = ROIs(((r-1)*80)+(1:80),:,:);
    UXX = UX.*qx.*mask.*ROI; UXY = UY.*qx.*mask.*ROI; UXZ = UZ.*qx.*mask.*ROI;
    UYX = UX.*qy.*mask.*ROI; UYY = UY.*qy.*mask.*ROI; UYZ = UZ.*qy.*mask.*ROI;
    UZX = UX.*qz.*mask.*ROI; UZY = UY.*qz.*mask.*ROI; UZZ = UZ.*qz.*mask.*ROI;

    UZX(UZX==0) = nan;UXX(UXX==0) = nan;UYX(UYX==0) = nan;
    UZZ(UZZ==0) = nan;UXZ(UXZ==0) = nan;UYZ(UYZ==0) = nan;
    UZY(UZY==0) = nan;UXY(UXY==0) = nan;UYY(UYY==0) = nan;

    figure;
    hold on
    quivtmp1 = quiver3(YY,XX,ZZ,UXX,UXY,UXZ,2,'r');
    quivtmp2 = quiver3(YY,XX,ZZ,UYX,UYY,UYZ,2,'color',[0 0.65 0]);
    quivtmp3 = quiver3(YY,XX,ZZ,UZX,UZY,UZZ,2,'b');
    hold off

%     figure;
%     hold on
%     im(tmpflip1(:,:,30))
%     quiver(XX(:,:,30),YY(:,:,30),UXX(:,:,30),UXY(:,:,30),3,'r'); axis equal;
%     quiver(XX(:,:,30),YY(:,:,30),UYX(:,:,30),UYY(:,:,30),3,'color',[0 .65 0]); axis equal;
%     quiver(XX(:,:,30),YY(:,:,30),UZX(:,:,30),UZY(:,:,S_z/2),3,'b'); axis equal;
%     hold off

    %XXp = permute(XX,[3 2 1]); ZZp1 = permute(ZZ,[3 2 1]); 
    %UXXp = permute(UXX,[3 2 1]); UXZp1 = permute(UXZ,[3 2 1]);
    %UYXp = permute(UYX,[3 2 1]); UYZp1 = permute(UYZ,[3 2 1]);
    %UZXp = permute(UZX,[3 2 1]); UZZp1 = permute(UZZ,[3 2 1]);
    
    %figure;
    %hold on
    %quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UXXp(:,:,S_y/2),UXZp1(:,:,S_y/2),3,'r'); axis equal;
    %quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UYXp(:,:,S_y/2),UYZp1(:,:,S_y/2),3,'color',[0 .65 0]); axis equal;
    %quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UZXp(:,:,S_y/2),UZZp1(:,:,S_y/2),3,'b'); axis equal;
    %hold off


    %YYp = permute(YY,[3 1 2]); ZZp2 = permute(ZZ,[3 1 2]);
    %UXYp = permute(UXY,[3 1 2]); UXZp2 = permute(UXZ,[3 1 2]);
    %UYYp = permute(UYY,[3 1 2]); UYZp2 = permute(UYZ,[3 1 2]);
    %UZYp = permute(UZY,[3 1 2]); UZZp2 = permute(UZZ,[3 1 2]);

    %figure;
    %hold on
    %quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UXYp(:,:,S_x/2+2),UXZp2(:,:,S_x/2+2),3,'r'); axis equal;
    %quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UYYp(:,:,S_x/2+2),UYZp2(:,:,S_x/2+2),3,'color',[0 .65 0]); axis equal;
    %quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UZYp(:,:,S_x/2+2),UZZp2(:,:,S_x/2+2),3,'b'); axis equal;
    %hold off
end

