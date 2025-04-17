function [Dr,Di]=tet_interp_disp(xyz,Ur,Ui,vox,mask,MagIm)
% interpolates the measured displacement field to the FE mesh
% Inputs: xyz = xyz positions of FE nodes (m)
%         Ur, Ui : Real and Imaginary displacement stacks
%         vox = size of each Ur and Ui voxel (m) [sx sy sz]
%         mask, MagIm (optional): Just to create a plot.


% pixel spacing
dx=vox(1);
dy=vox(2);
dz=vox(3);

% displacements - stack without top&bottom slices
[nx,ny,nslice,nd]=size(Ur);
clear A P;

% 'pad' displ
pDr=zeros(nx,ny,nslice+2,nd);
pDi=zeros(nx,ny,nslice+2,nd);

pDr(:,:,1,:)=Ur(:,:,1,:);
pDr(:,:,2:nslice+1,:)=Ur;
pDr(:,:,nslice+2,:)=Ur(:,:,end,:);

pDi(:,:,1,:)=Ui(:,:,1,:);
pDi(:,:,2:nslice+1,:)=Ui;
pDi(:,:,nslice+2,:)=Ui(:,:,end,:);
%% Create mesh grid to interpolate displacements and display

% z dir includes the 'extra' slices in the nslice dimension
[xx,yy,zz]=meshgrid(0:dx:(nx-1)*dx,0:dy:(ny-1)*dy,-dz:dz:(nslice)*dz);

MINx = min(xx(:));MAXx = max(xx(:));
MINy = min(yy(:));MAXy = max(yy(:));
MINz = min(zz(:));MAXz = max(zz(:));

dim = [MINx MINy MINz; ...
       MINx MAXy MINz; ...
       MAXx MINy MINz; ...
       MAXx MAXy MINz; ...
       MINx MINy MAXz; ...
       MINx MAXy MAXz; ...
       MAXx MINy MAXz; ...
       MAXx MAXy MAXz];



%% Plot MRE First & Last Images

if (nargin==6)
    plot3(dim(:,1),dim(:,2),dim(:,3),'.r');
    grid on;
    hold on;
    
    MagIm = MagIm.*mask;

    X=[1,nx;1,nx]*dx;
    Y=[1,1;ny,ny]*dy;
    Z=[1,1;1,1]*dz;
    surface(X,Y,Z,'CData',double(MagIm(:,:,1)'),'FaceColor','texturemap');
    colormap(gray);drawnow;

    Z=[1,1;1,1]*dz*(nslice);
    surface(X,Y,Z,'CData',double(MagIm(:,:,end)'),'FaceColor','texturemap');
    colormap(gray);drawnow;

    %% Display FE mesh
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'g.')
%     bel = load([series_name,'.bel']);
%     trisurf(bel(:,2:4),xyz(:,2),xyz(:,3),xyz(:,4),'facecolor','g','EdgeColor','none');
%     axis equal;xlabel('\bfX'),ylabel('\bfY'),zlabel('\bfZ');
end
%saveas(gcf,'mesh-validation.jpg','jpeg');%close;

%% Interpolate displacements from grid to the FE mesh
Dr(:,1)=interp3(xx,yy,zz,permute(pDr(:,:,:,1),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');
Dr(:,2)=interp3(xx,yy,zz,permute(pDr(:,:,:,2),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');
Dr(:,3)=interp3(xx,yy,zz,permute(pDr(:,:,:,3),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');

Di(:,1)=interp3(xx,yy,zz,permute(pDi(:,:,:,1),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');
Di(:,2)=interp3(xx,yy,zz,permute(pDi(:,:,:,2),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');
Di(:,3)=interp3(xx,yy,zz,permute(pDi(:,:,:,3),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'*cubic');

end
