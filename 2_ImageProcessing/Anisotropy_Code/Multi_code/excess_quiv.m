%extra quiver code

%zeromask = ones(size(XX));
%for i = 1:S_x
%    for j = 1:S_y
%        for k = 1:S_z
%            if UXX1(i,j,k) == 0 || UYX1(i,j,k) == 0 || UZX1(i,j,k) == 0 || UXY1(i,j,k) == 0 || UXZ1(i,j,k) == 0 || UYY1(i,j,k) == 0 || UYZ1(i,j,k) == 0 || UZZ1(i,j,k) == 0 || UZY1(i,j,k) == 0
%                zeromask(i,j,k) = 0;
%            elseif isnan(UXX1(i,j,k)) == 0 || isnan(UYX1(i,j,k)) == 0 || isnan(UZX1(i,j,k)) == 0 || isnan(UXY1(i,j,k)) == 0 || isnan(UXZ1(i,j,k)) == 0 || isnan(UYY1(i,j,k)) == 0 || isnan(UYZ1(i,j,k)) == 0 || isnan(UZZ(i,j,k)) == 0 || isnan(UZY(i,j,k)) == 0
%                zeromask(i,j,k) = 0;
%            end
%        end
%    end
%end

%UX1 = Output(:,:,:,1); UX1 = UX1(1:nskip:xy,1:nskip:xy,1:nskip:z); 
%UY1 = Output(:,:,:,2); UY1 = UY1(1:nskip:xy,1:nskip:xy,1:nskip:z);  
%UZ1 = Output(:,:,:,3); UZ1 = UZ1(1:nskip:xy,1:nskip:xy,1:nskip:z); 

%UXX1 = UX1.*qx.*mask; UXY1 = UY1.*qx.*mask; UXZ1 = UZ1.*qx.*mask;
%UYX1 = UX1.*qy.*mask; UYY1 = UY1.*qy.*mask; UYZ1 = UZ1.*qy.*mask;
%UZX1 = UX1.*qz.*mask; UZY1 = UY1.*qz.*mask; UZZ1 = UZ1.*qz.*mask;
%cUXX = UXX1.*zeromask; cUXY = UXY1.*zeromask; cUXZ = UXZ1.*zeromask;
%cUYX = UXX1.*zeromask; cUYY = UXY1.*zeromask; cUYZ = UXZ1.*zeromask;
%cUZX = UXX1.*zeromask; cUZY = UZY1.*zeromask; cUZZ = UXZ1.*zeromask;

%figure;
%hold on
%for i = 25:45
%    for j = 25:45
%        for k = 20:30
%            [i,j,k]
%            quiver3(XX(i,j,k),YY(i,j,k),ZZ(i,j,k),UXX(i,j,k),UXY(i,j,k),UXZ(i,j,k),2,'color',[cUXX(i,j,k) cUXY(i,j,k) cUXZ(i,j,k)]);
%            quiver3(XX(i,j,k),YY(i,j,k),ZZ(i,j,k),UYX(i,j,k),UYY(i,j,k),UYZ(i,j,k),2,'color',[cUYX(i,j,k) cUYY(i,j,k) cUYZ(i,j,k)]);
%            quiver3(XX(i,j,k),YY(i,j,k),ZZ(i,j,k),UZX(i,j,k),UZY(i,j,k),UZZ(i,j,k),2,'color',[cUZX(i,j,k) cUZY(i,j,k) cUZZ(i,j,k)]);
%        end
%    end
%end
%hold off

