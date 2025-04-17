function [regout]=tet_nearest_interp(xyz,stack,vox)

% Pad stack
[nx,ny,nslice]=size(stack);

pstack=zeros(nx,ny,nslice+2);

pstack(:,:,1)=stack(:,:,1);
pstack(:,:,2:nslice+1)=stack;
pstack(:,:,nslice+2,:)=stack(:,:,end);

[xx,yy,zz]=meshgrid(0:vox(1):(nx-1)*vox(1),0:vox(2):(ny-1)*vox(2),-vox(3):vox(3):(nslice)*vox(3));

% xx=permute(xx,[2 1 3]);
% yy=permute(yy,[2 1 3]);


% Permute input stack to fit the (:,1,1)=x, (1,:,1)=y assumption
regout(:,1)=interp3(xx,yy,zz,permute(pstack(:,:,:,1),[2 1 3]),xyz(:,1),xyz(:,2),xyz(:,3),'nearest');

end



