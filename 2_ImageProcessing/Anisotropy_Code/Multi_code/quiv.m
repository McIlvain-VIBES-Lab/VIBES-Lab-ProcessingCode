function quiv(Output,xy,z,a,scaling)
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

for i = 3:S_x-2
    for j = 3:S_y-2
        for k = 3:S_z-2
            if UX(i+2,j,k) == 0 || UX(i-2,j,k) == 0 || UY(i+2,j,k) == 0 || UY(i-2,j,k) == 0 || UZ(i+2,j,k) == 0 || UZ(i-2,j,k) == 0 
                zerom(i,j,k) = 1;
            elseif UX(i,j+2,k) == 0 || UX(i,j-2,k) == 0 || UY(i,j+2,k) == 0 || UY(i,j-2,k) == 0 || UZ(i,j+2,k) == 0 || UZ(i,j-2,k) == 0 
                zerom(i,j,k) = 1;
            elseif UX(i,j,k+2) == 0 || UX(i,j,k-2) == 0 || UY(i,j,k+2) == 0 || UY(i,j,k-2) == 0 || UZ(i,j,k+2) == 0 || UZ(i,j,k+2) == 0 
                zerom(i,j,k) = 1;
            end
        end
    end
end

oneslice = ones(S_x,S_y); mask = onem - zerom;
mask(:,:,1) = mask(:,:,1) - oneslice; mask(:,:,2) = mask(:,:,2) - oneslice;
mask(:,:,S_z) = mask(:,:,S_z) - oneslice; mask(:,:,S_z-1) = mask(:,:,S_z-1) - oneslice;

XX = X(1:nskip:xy,1:nskip:xy,1:nskip:z);
YY = Y(1:nskip:xy,1:nskip:xy,1:nskip:z);
ZZ = Z(1:nskip:xy,1:nskip:xy,1:nskip:z);

UXX = UX.*qx.*mask; UXY = UY.*qx.*mask; UXZ = UZ.*qx.*mask;
UYX = UX.*qy.*mask; UYY = UY.*qy.*mask; UYZ = UZ.*qy.*mask;
UZX = UX.*qz.*mask; UZY = UY.*qz.*mask; UZZ = UZ.*qz.*mask;

%for i = 1:S_x
%   for j = 1:S_y
%        for k = 1:S_z
%            if UX(i,j,k) == 0 || UY(i,j,k) == 0 || UZ(i,j,k) == 0
%               UX(i,j,k) = 0;
%               UY(i,j,k) = 0;
%               UZ(i,j,k) = 0;
%            end
%        end
%   end
%end
zeromask = ones(size(XX));
for i = 1:S_x
    for j = 1:S_y
        for k = 1:S_z
            if UXX(i,j,k) == 0 || UYX(i,j,k) == 0 || UZX(i,j,k) == 0 || UXY(i,j,k) == 0 || UXZ(i,j,k) == 0 || UYY(i,j,k) == 0 || UYZ(i,j,k) == 0 || UZZ(i,j,k) == 0 || UZY(i,j,k) == 0
                zeromask(i,j,k) = 0;
            elseif isnan(UXX1(i,j,k)) == 0 || isnan(UYX1(i,j,k)) == 0 || isnan(UZX1(i,j,k)) == 0 || isnan(UXY1(i,j,k)) == 0 || isnan(UXZ1(i,j,k)) == 0 || isnan(UYY1(i,j,k)) == 0 || isnan(UYZ1(i,j,k)) == 0 || isnan(UZZ(i,j,k)) == 0 || isnan(UZY(i,j,k)) == 0
                zeromask(i,j,k) = 0;
            end
        end
    end
end

UZX(UZX==0) = nan;UXX(UXX==0) = nan;UYX(UYX==0) = nan;
UZZ(UZZ==0) = nan;UXZ(UXZ==0) = nan;UYZ(UYZ==0) = nan;
UZY(UZY==0) = nan;UXY(UXY==0) = nan;UYY(UYY==0) = nan;


%size(ZZ)
%size(UX)

%figure;quivtmp1 = quiver3(XX,YY,ZZ,UXX,UXY,UXZ,2,'k');
%figure;quivtmp2 = quiver3(XX,YY,ZZ,UYX,UYY,UYZ,2,'b');
%figure;quivtmp3 = quiver3(XX,YY,ZZ,UZX,UZY,UZZ,2,'r');

figure;
hold on
quiver3(XX,YY,ZZ,UXX,UXY,UXZ,2,'r'); axis equal;
quiver3(XX,YY,ZZ,UYX,UYY,UYZ,2,'color',[0 .65 0]); axis equal;
quiver3(XX,YY,ZZ,UZX,UZY,UZZ,2,'b'); axis equal;
hold off

%figure;quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UXX(:,:,S_z/2),UXY(:,:,S_z/2),2,'k');
%figure;quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UYX(:,:,S_z/2),UYY(:,:,S_z/2),2,'b');
%figure;quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UZX(:,:,S_z/2),UZY(:,:,S_z/2),2,'r');

figure;
hold on
quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UXX(:,:,S_z/2),UXY(:,:,S_z/2),3,'r'); axis equal;
quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UYX(:,:,S_z/2),UYY(:,:,S_z/2),3,'color',[0 .65 0]); axis equal;
quiver(XX(:,:,S_z/2),YY(:,:,S_z/2),UZX(:,:,S_z/2),UZY(:,:,S_z/2),3,'b'); axis equal;
hold off

XXp = permute(XX,[3 2 1]); ZZp1 = permute(ZZ,[3 2 1]); 
UXXp = permute(UXX,[3 2 1]); UXZp1 = permute(UXZ,[3 2 1]);
UYXp = permute(UYX,[3 2 1]); UYZp1 = permute(UYZ,[3 2 1]);
UZXp = permute(UZX,[3 2 1]); UZZp1 = permute(UZZ,[3 2 1]);
figure;
hold on
quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UXXp(:,:,S_y/2),UXZp1(:,:,S_y/2),3,'r'); axis equal;
quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UYXp(:,:,S_y/2),UYZp1(:,:,S_y/2),3,'color',[0 .65 0]); axis equal;
quiver(XXp(:,:,S_y/2),ZZp1(:,:,S_y/2),UZXp(:,:,S_y/2),UZZp1(:,:,S_y/2),3,'b'); axis equal;
hold off


YYp = permute(YY,[3 1 2]); ZZp2 = permute(ZZ,[3 1 2]);
UXYp = permute(UXY,[3 1 2]); UXZp2 = permute(UXZ,[3 1 2]);
UYYp = permute(UYY,[3 1 2]); UYZp2 = permute(UYZ,[3 1 2]);
UZYp = permute(UZY,[3 1 2]); UZZp2 = permute(UZZ,[3 1 2]);

figure;
hold on
quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UXYp(:,:,S_x/2+2),UXZp2(:,:,S_x/2+2),3,'r'); axis equal;
quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UYYp(:,:,S_x/2+2),UYZp2(:,:,S_x/2+2),3,'color',[0 .65 0]); axis equal;
quiver(YYp(:,:,S_x/2+2),ZZp2(:,:,S_x/2+2),UZYp(:,:,S_x/2+2),UZZp2(:,:,S_x/2+2),3,'b'); axis equal;
hold off

