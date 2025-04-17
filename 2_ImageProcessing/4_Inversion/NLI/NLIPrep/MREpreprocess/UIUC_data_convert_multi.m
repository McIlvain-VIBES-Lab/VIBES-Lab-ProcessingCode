% Convert Curtis johnsons MRE data (unversity of Illinois Urbana Champaign)
% to our format (old MRE_3DMotionData format, which is then converted to the new
% NLI_MRE_Input format.

% Our format:
% File MRE-3DMotionData.mat contains
%   A : real valued [size nx ny nz 3] %Motion Amplitude (4th dimension is (u,v,w))
%   P : real valued [size nx ny nz 3] %Phase angle of complex motions, so
%                                      that Ureal=A.*cos(P) and Uimag=A.*sin(P)
%   MagIm : real valued [size nx ny nz 3] % MR magnitude image or some 
%                                           other anatomical image
%
% File HeaderData.mat contains
%   DirIndex % This is a confusing variable. It contains the required
%              transformation matrices from the matlab indices to the xyz
%              coordinate system used in the mesh. It gets more confusing:
%              DirIndex(1:3,1:3) transforms the physical directions, and
%              DirIndex(1:3,4:6) transforms the motion directions. I think 
%              these two are redundant, i.e. you only need one. Although I 
%              am not sure, i find all this direction stuff horribly 
%              confusing. We have always used DirIndex(1:3,1:3)=eye(3),
%              which implies that MagIm(:,1,1) is a line in the x
%              direction, MagIm(1,:,1) is the y direction, and MagIm(1,1,:)
%              is the z direction. Then, as an example, if we determined 
%              a positive number in Ureal(:,:,:,1) is a motion in the +ve y
%              direction, a positive number in Ureal(:,:,:,2) is a motion 
%              in the -ve x direction, and  positive number in 
%              Ureal(:,:,:,3) is a motion in the +ve z direction, then 
%              Dirindex(1:3,4:6) = [0 1 0;-1 0 0;0 0 1]. This can be used 
%              to transform the displacmeents like this:
%              % Transform to Nelx3 array for transformation
%              Uvec123(:,1)=reshape(Ur(:,:,:,1),[nx*ny*nz,1]);
%              Uvec123(:,2)=reshape(Ur(:,:,:,2),[nx*ny*nz,1]);
%              Uvec123(:,3)=reshape(Ur(:,:,:,3),[nx*ny*nz,1]);
%              Uvecxyz=Uvec123*DirIndex(1:3,4:6)
%              % Transform back to image stack 
%              Ur_xyz(:,1)=reshape(Uvecxyz(:,:,:,1),[nx ny nz]);
%              Ur_xyz(:,2)=reshape(Uvecxyz(:,:,:,2),[nx ny nz]);
%              Ur_xyz(:,3)=reshape(Uvecxyz(:,:,:,3),[nx ny nz]);
%              Then, DirIndex(4,1:3) is the voxel size in each of the 3
%              directions, in mm. DirIndex(4,4:6) is not used.
%  freqHz % Frequnecy in Hz
%
% File Mask.mat contains
%   mask : size [nx ny nz]
%   
clear all
%fname=getfilename('Data file(s) to convert >>','*.mat');

d=dir('*.mat');
for ii=1:length(d)
    disp(['oo    File ' int2str(ii) ' :: ' d(ii).name])
end
disp('oo    Select multiple disp files for multi-set inversion e.g. [4 5]')
fnums=input('oo    Enter Number of file(s) to convert (default=[1]) >> ');

load(d(fnums(1)).name);
Ur=zeros([size(Xmotion) 3 length(fnums)]);
Ui=zeros([size(Xmotion) 3 length(fnums)]);

for ii=1:length(fnums)
    
    % Tidy up inconsistent mask names in supplied files (should be mask,
    % sometimes new_mask, sometimes eroded_mask
    clear new_mask eroded_mask
    inputfiles(ii).name=d(fnums(ii)).name;
    load(inputfiles(ii).name)
    if(exist('new_mask','var'))        
        mask=new_mask;
    elseif(exist('eroded_mask','var'))
        mask=eroded_mask;
    end
    
    % Check all masks are the same. Use first set if not. 
    if(ii==1)
        mask1=mask;
        AnatomicalMRI=t2stack; % Save the first T2 as MagIm
    else
        if(sum(abs(mask(:)-mask1(:)))~=0)
            disp('WARNING!!!! Masks appear to be different, must be the same for multi-set inversion')
            disp('Setting mask to the intersetion of all masks')
            mask=mask&mask1;
        end
    end
    
    % Hard coded transoframtions, assuming axial, HFS, Phase Encode =
    % whatever it is for curtis's standard brain imaging. Saggital,
    % coronal, FFS, different phase encode directions will probably make
    % this wrong.
    rowswap=[0 1 0; 1 0 0 ; 0 0 1];
    Motion2Image=rowswap*[-1 0 0;0 1 0;0 0 1];
    RHcoord=true;

    voxsize_mm(1)=mreParams.FOVx/mreParams.nx;
    voxsize_mm(2)=mreParams.FOVy/mreParams.ny;
    voxsize_mm(3)=mreParams.FOVz/mreParams.nz;

    
    freqHz(ii)=mreParams.freq;
    
    

    Ur(:,:,:,1,ii)=real(Xmotion);
    Ui(:,:,:,1,ii)=imag(Xmotion);
    Ur(:,:,:,2,ii)=real(Ymotion);
    Ui(:,:,:,2,ii)=imag(Ymotion);
    Ur(:,:,:,3,ii)=real(Zmotion);
    Ui(:,:,:,3,ii)=imag(Zmotion);
end
dname=[inputfiles(1).name(1:end-4) '_' int2str(length(fnums)) 'set'];
sv=false;
if(~exist(dname,'dir'))
    mkdir(dname)
    sv=true;
else
    disp(['WARNING:: directory ' dname ' Exists'])
    qn=input('overwrite files <y/n> ??','s');
    if(strcmp(qn,'y'))
        sv=true;
    end
end
if(sv)
    save([dname '/NLI_MRE_Input.mat'],'Ur','Ui','voxsize_mm','RHcoord','Motion2Image','AnatomicalMRI','freqHz','inputfiles')
    save([dname '/Mask.mat'],'mask')    
    % Process all the different SPR name files we have got 
    % Can only have one SPR regs file. Check to see if its a variable. 
    load(inputfiles(1).name);
    if(exist('SPR','var'))
        sspr=size(SPR);
        regs=zeros(sspr(1:3));
        for jj=1:sspr(4)
            regs(SPR(:,:,:,jj)==1)=jj;
        end
        save([dname '/SPRregs.mat'],'regs')
    end
    % Check to see if it is in a supplied file from UIUC
    if(exist('SPR.mat','file'))
        load SPR.mat
        regs=zeros([size(SPR,1) size(SPR,2) size(SPR,3)]);
        for kk=1:size(SPR,4)
            regs(SPR(:,:,:,kk)==1)=kk;
        end
        save([dname '/SPRregs.mat'],'regs')
    end
    
    % Check to see if we have made one with any of my codes (eg
    % Stack_Segmentation.m, 
    f=[dir('*SPRregs.mat')];
    if(length(f)>0)
        copyfile(f(1).name,dname);
    end
    
    if(exist('segmentation.mat','file'))
        copyfile('segmentation.mat',dname);
    end
    
    % Check for DTI file
    if(exist('DTI.mat','file'))
        copyfile('DTI.mat',dname);
    end

end

