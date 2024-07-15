function convdata(voxres,freq,wimg,eimg,t2stack,mask)

% FROM: Matthew D. McGarry, Dartmouth College
% MODIFIED by Curtis L. Johnson for use in MRE codes



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



% clear all
% fname=getfilename('Data file to convert >>','*.mat');

% load(fname)
DirIndex=zeros(4,6);

DirIndex(1:3,1:3)=eye(3);
DirIndex(1:3,4:6)=eye(3);
DirIndex(4,1)=voxres(1);
DirIndex(4,2)=voxres(2);
DirIndex(4,3)=voxres(3);

MagIm=t2stack;
freqHz=freq;

A(:,:,:,1)=abs(wimg(:,:,:,3));
P(:,:,:,1)=angle(wimg(:,:,:,3));
A(:,:,:,2)=abs(wimg(:,:,:,2));
P(:,:,:,2)=angle(wimg(:,:,:,2));
A(:,:,:,3)=abs(wimg(:,:,:,1));
P(:,:,:,3)=angle(wimg(:,:,:,1));

% switching the order of directions for errormap
ErrorMap(:,:,:,1) = eimg(:,:,:,3);
ErrorMap(:,:,:,2) = eimg(:,:,:,2);
ErrorMap(:,:,:,3) = eimg(:,:,:,1);
stdform=true;

% dname=fname(1:end-4);
dname = 'dcfiles';
sv=true;
if(~exist(dname,'dir'))
    mkdir(dname)
    % sv=true;
% else
%     disp(['WARNING:: directory ' dname ' Exists'])
%     qn=input('overwrite files <y/n> ??','s');
%     if(strcmp(qn,'y'))
%         sv=true;
%     end
end
if(sv)
    save([dname '/HeaderData.mat'],'DirIndex','freqHz')
    save([dname '/MRE_3DMotionData.mat'],'A','P','MagIm')
    save([dname '/Mask.mat'],'mask')
    save([dname '/ErrorMap.mat'],'ErrorMap','stdform')
end

% You would need to add something here to generate ErrorMap.mat to make my
% StrainSNR.m code work.

