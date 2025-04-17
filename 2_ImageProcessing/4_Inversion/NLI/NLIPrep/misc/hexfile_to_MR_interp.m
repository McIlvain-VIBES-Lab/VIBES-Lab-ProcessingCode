% Interpolates data from a hexahedral mesh onto the voxels of an MR stack.


function [stack]=hexfile_to_MR_interp(hexf,nodf,dim,res)
% hexf is a file that contains a variable on hexahedral nodes to be
% interpolated to an MR stack, where dim=size(stack) and res = the
% voxel size of the stack.
% Inputs: hexf: file to interpolate
%         nodf: Appropriate node file
%         If dim and res are not given, fuction will prompt for a MRE file
%         to take them from. Or, supply them:
%         dim: the size of the MR stack [nx ny nz]
%         res: Voxel size in m. [rx ry rz]
% This function assumes the same conventions as the MRE preprocessing for
% the origin location and x/y/z axis definitions.
% It also assumes that the nodes are in a rgular grid (which they will be
% for any hex mesh in the current MRE pipeline.


if(nargin<1)
    [hexf,hexpath]=uigetfile('*','Select hexahedral variable file to interpolate');
    hexval=load(fullfile(hexpath,hexf));
else
    hexval=load(hexf);
end

if(nargin<2)
    [nodf,nodpath]=uigetfile('*.nod','Select appropriate hexahedral node file');
    nod=load(fullfile(nodpath,nodf));
else
    nod=load(nodf);
end

if(size(nod,1)~=size(hexval,1))
    error('hexf and nodf must have same number of rows')
end

if(nargin<3)
    [MRf,MRpath]=uigetfile('*','Select MRE file to interpolate back to');
    load(fullfile(MRpath,MRf));
    if(exist('A','var')) % Old MRE_3DMotionData format
        dim=[size(A,1) size(A,2) size(A,3)];
        load(fullfile(MRpath,'HeaderData.mat'));
        res=DirIndex(4,1:3)/1000;
    elseif(exist('Ur','var'))
        dim=[size(Ur,1) size(Ur,2) size(Ur,3)];
        res=voxsize_mm/1000;
        
    end
end

% Strip node numbers
hexval=hexval(:,2:end);
x=nod(:,2);
y=nod(:,3);
z=nod(:,4);
sense=nod(:,5);
nset=size(hexval,2);
nn=size(nod,1);

stack=zeros([dim nset]);

% Sort the hex nodes into a grid;
unqx=sort(unique(x));
unqy=sort(unique(y));
unqz=sort(unique(z));
reshex=[mean(diff(unqx)) mean(diff(unqy)) mean(diff(unqz))];

indx=round((x-min(x))/reshex(1))+1;
indy=round((y-min(y))/reshex(1))+1;
indz=round((z-min(z))/reshex(1))+1;

% indx=indx-min(indx)+1;
% indy=indy-min(indy)+1;
% indz=indz-min(indz)+1;

hexstack=nan(max(indx),max(indy),max(indz),nset);

% Same coordinate assumptions as MRE preprocessing code
Xmr=0:res(1):res(1)*(dim(1)-1);
Ymr=0:res(2):res(2)*(dim(2)-1);
Zmr=0:res(3):res(3)*(dim(3)-1);

sumsense=sum(sense); % Check that the sensitivities aresnt all zero (pre MREv9p33)
for ii=1:nset
    for jj=1:nn
        if(sense(jj)||(sumsense==0)) % Leave as nan if node outside the disp mesn (v9.33 and above). 
            hexstack(indx(jj),indy(jj),indz(jj),ii)=hexval(jj,ii);
        end
    end
    stack(:,:,:,ii)=LinearInterp3D_withnans(unqx,unqy,unqz,hexstack(:,:,:,ii),Xmr,Ymr,Zmr);
end

end