function UIUC_data_convert_mcilvain(SubjectName)
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
%clear all
%fname=getfilename('Data file(s) to convert >>','*.mat');

d=(sprintf('%s.mat',SubjectName));
% for ii=1:length(d);
%     disp(['oo    File ' int2str(ii) ' :: ' d(ii).name])
% end
fnums=1; %input('oo    Enter Number of file(s) to convert (default=[1]) >> ');


for ii=1:length(fnums)
    clear new_mask
    fname=d;
    load(fname)
    if(exist('new_mask','var'))
        mask=new_mask;
    elseif(exist('eroded_mask','var'))
        mask=eroded_mask;
    end
    
    DirIndex=zeros(4,6);

    DirIndex(1:3,1:3)=eye(3);
    %DirIndex(1:3,4:6)=[0 1 0;-1 0 0;0 0 1];
    DirIndex(1:3,4:6)=[-1 0 0;0 1 0;0 0 1];

    DirIndex(4,1)=mreParams.FOVx/mreParams.nx;
    DirIndex(4,2)=mreParams.FOVy/mreParams.ny;
    DirIndex(4,3)=mreParams.FOVz/mreParams.nz;

    MagIm=t2stack;
    freqHz=mreParams.freq;


    scalfac=10;
    Xmotion=Xmotion*scalfac;
    Ymotion=Ymotion*scalfac;
    Zmotion=Zmotion*scalfac;
    
    clear A P

    A(:,:,:,1)=abs(Xmotion);
    P(:,:,:,1)=angle(Xmotion);
    A(:,:,:,2)=abs(Ymotion);
    P(:,:,:,2)=angle(Ymotion);
    A(:,:,:,3)=abs(Zmotion);
    P(:,:,:,3)=angle(Zmotion);

    dname=fname(1:end-4);
 %   sv=false;
 %   if(~exist(dname,'dir'))
        mkdir(dname)
        sv=true;
 %   else
 %       disp(['WARNING:: directory ' fname ' Exists'])
 %       qn=input('overwrite files <y/n> ??','s');
 %       if(strcmp(qn,'y'))
 %           sv=true;
 %       end
 %   end
    if(sv)
        save([dname '/HeaderData.mat'],'DirIndex','freqHz','scalfac')
        save([dname '/MRE_3DMotionData.mat'],'A','P','MagIm')
        save([dname '/Mask.mat'],'mask')
        
        if(exist('SPR','var'))
            sspr=size(SPR);
            regs=zeros(sspr(1:3));
            for jj=1:sspr(4);
                regs(SPR(:,:,:,jj)==1)=jj;
            end
            save([dname '/SPRregs.mat'],'regs')
        end
        if(exist('segmentation.mat','file'))
            copyfile('segmentation.mat',dname);
        end
        if(exist('SPR.mat','file'))
            load SPR.mat
            regs=zeros([size(SPR,1) size(SPR,2) size(SPR,3)]);
            for kk=1:size(SPR,4)
                regs(SPR(:,:,:,kk)==1)=kk;
            end
            save([dname '/SPRregs.mat'],'regs')
        end
		f=[dir('*SPRregs.mat')];
		if(length(f)>0)
		    copyfile(f(1).name,dname);
        end
        
        cd(dname);
        MRE_MotionData_convert;
        delete('MRE_3DMotionData.mat');
        delete('HeaderData.mat');
        cd ..
        
    end
end
% You would need to add something here to generate ErrorMap.mat to make my
% StrainSNR.m code work.

end