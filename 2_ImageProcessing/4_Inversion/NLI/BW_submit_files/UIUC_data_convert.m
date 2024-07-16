% Convert Curtis johnsons MRE data (unversity of Illinois Urbana Champaign)
% to our format

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
% File ErrorMap.mat contains
%   ErrorMap : real valued [size nx ny nz 3]
%     An estimate of the displacement error in each of the 3 directions. The 
%     simplest way to generate it is to use the standard devaition of the
%     misfits between the fitted amplitude and the data at each phase offset.
%     For example, if rawdata was a [nx ny nz 3 nph] array containing the
%     nph raw phase images in each of the 3 directions. The fundamental
%     mode of the y directed motion of voxel [10 10 10] would be given by
%     F=fft(rawdata(10,10,10,2,:)), Ucmplx(10,10,10,2)=2/nph*F(2). Then, you 
%     could estimate the uncertianty of the amplitude of that voxel by: 
%     t= 0:2*pi/nph:(2*nph-1)/nph*pi;
%     ErrorMap(10,10,10,2) = std(real(2/nph*F(2)*exp(1i*t)) - rawdata(10,10,10,2,:))
%   stdform : set to true if using the form of ErrorMap described above.
%     This is an indicator I added when I moved from using a errormap defined 
%     by mean(misfit) to one defined by std(misfit). I apply a conversion factor 
%     to the old mean(misfit) form so I didnt have to repeat all the motion 
%     unwrapping etc for our old data. mean(abs(r)) ~= 0.7979*std(r) for
%     normally distributed r. 
%
%   There is another, possibly better way of calcualting strain noise,
%   using propagation of uncertianties to map the noise in the raw MR image
%   into noise in the octahedral shear strain. I have used propogation of
%   uncertainties to go from noise in the raw MR images to noise in the
%   motion, but I am too scared to try the octahedral shear equation and
%   all it's numerical derivatives. I will attach the unfinshed section of 
%   my thesis proposal with the derivation. This stuff is really pushing 
%   the limits of my understanding of MRI, so if you happen to read it and
%   see anything wrong let me know. I sometimes use this measure as a 
%   ErrorMap, I call the variable ErrorMapPropofErrs, you can see where it 
%   looks for it in the StrainSNR.m code. It ends up as 
%   ErrorMapPropofErrs(i,j,k)=MRnoise*sqrt(2/nph)*abs(MRmagnitude(i,j,k))
%   
%   



clear all
%fname=getfilename('Data file(s) to convert >>','*.mat');

d=dir('*.mat');
for ii=1:length(d);
    disp(['oo    File ' int2str(ii) ' :: ' d(ii).name])
end
fnums=input('oo    Enter Number of file(s) to convert (default=[1]) >> ');


for ii=1:length(fnums)

    fname=d(fnums(ii)).name;
    load(fname)
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


    A(:,:,:,1)=abs(Xmotion);
    P(:,:,:,1)=angle(Xmotion);
    A(:,:,:,2)=abs(Ymotion);
    P(:,:,:,2)=angle(Ymotion);
    A(:,:,:,3)=abs(Zmotion);
    P(:,:,:,3)=angle(Zmotion);

    dname=fname(1:end-4);
    sv=false;
    if(~exist(dname,'dir'))
        mkdir(dname)
        sv=true;
    else
        disp(['WARNING:: directory ' fname ' Exists'])
        qn=input('overwrite files <y/n> ??','s');
        if(strcmp(qn,'y'))
            sv=true;
        end
    end
    if(sv)
        save([dname '/HeaderData.mat'],'DirIndex','freqHz','scalfac')
        save([dname '/MRE_3DMotionData.mat'],'A','P','MagIm')
        save([dname '/Mask.mat'],'mask')
        
        if(exist('SPR','var'))
            sspr=size(SPR);
            regs=zeros(sspr(1:3));
	    if(length(sspr)==4)
	        for jj=1:sspr(4);
                    regs(SPR(:,:,:,jj)==1)=jj;
                end
	    else
		regs(SPR(:,:,:)==1)=1;
	    end
            save([dname '/SPRregs.mat'],'regs')
        end          
        
    end
end
% You would need to add something here to generate ErrorMap.mat to make my
% StrainSNR.m code work.

