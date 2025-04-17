% Convert WashU MRE data (unversity of Illinois Urbana Champaign)
% to our format (The new NLI_MRE_Input format.)

% Our format:
% File MRE-3DMotionData.mat contains
%   Ur : real valued [size nx ny nz 3] % Realmotion amplitude (4th dimension is (u,v,w))
%   Ui : real valued [size nx ny nz 3] % Phase angle of complex motions, so
%                                        that Ureal=A.*cos(P) and Uimag=A.*sin(P)
%   AnatomicalMRI : real valued [size nx ny nz] % MR magnitude image or some 
%                                           other anatomical image
%
% File HeaderData.mat contains
%   Motion2Image % This is a confusing variable. It contains the required
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
[n1,n2,n3,nph,ndir]=size(PCWU_BMR);
Ur=zeros([n1 n2 n3 3 length(fnums)]);
Ui=zeros([n1 n2 n3 3 length(fnums)]);

if(exist('Freq','var'))
    ActFreq=Freq;
end
if(~exist('ActFreq','var'))
    disp('Important: ActFreq not present, Frequency currently hardcoded to')
    HardCodedFrequency=100
    ActFreq=HardCodedFrequency;
end

for ii=1:length(fnums)
    
    inputfiles(ii).name=d(fnums(ii)).name;
    % Check all masks are the same. Use first set if not. 
    if(ii==1)
        mask1=mask;
        AnatomicalMRI=MAGmean; % Save the first anatomical MRI
    else
        if(sum(abs(mask(:)-mask1(:)))~=0)
            disp('WARNING!!!! Masks appear to be different, must be the same for multi-set inversion')
            disp('Setting mask to the intersetion of all masks')
            mask=mask&mask1;
        end
    end
    
    % Hard coded transformations, From Ruth:
%     (:,:,:,1:8,1) refers to displacements (encoded as phase in radians)  in the image read-out direction, with positive displacement in the positive read-out direction
%     (:,:,:,1:8,2) refers to  displacements in the image phase-encode direction, with positive displacement in the direction of positive phase-encode direction
%     (:,:,:,1:8,3) refers to displacement in the image slice direction with displacement in the direction of increasing slice number being positive.
%     The fourth index is increasing temporal phase with intervals of pi/4.
% 
%     With respect to the spatial dimensions,
%     (1,1,1) would be the upper left corner of the first image in the stack.
%     Traversing across a row, i.e. (1,1:120,1), would be increasing values in the image read-out direction.
%     Traversing down a column i.e. 1:120,1,1) would be increasing values in the image phase-encode direction
%     Traversing along the slices i.e. (:,:,1:48) would be increasing values in the slice direction.
% 
%     We use the same system for phantom, minipig and human data.
    
    
    Motion2Image=eye(3);
    RHcoord=true; % Set this to false for image orientations that do not produce a RH coordinate system. In my experience it doesnt change the inversion, but eliminates a potential source of error 
    
    if(length(voxsize)==3)
        voxsize_mm=voxsize;
    else % Assume isotropic voxels
        voxsize_mm(1)=voxsize;
        voxsize_mm(2)=voxsize;
        voxsize_mm(3)=voxsize;
    end
    
    freqHz(ii)=ActFreq;
    FFTmotion=2/nph*PC2micron*fft(PCWU_BMR,[],4);   % Should give amplitudes in units of microns. Gets normalized to average amilitude of 1mm in NLI so it doesnt matter provided we are not using nonlinear models.

    Ur(:,:,:,1,ii)=real(FFTmotion(:,:,:,2,2)); % Column, phase encode, 2 
    Ui(:,:,:,1,ii)=imag(FFTmotion(:,:,:,2,2));
    Ur(:,:,:,2,ii)=real(FFTmotion(:,:,:,2,1)); % row, readout, 1
    Ui(:,:,:,2,ii)=imag(FFTmotion(:,:,:,2,1));
    Ur(:,:,:,3,ii)=real(FFTmotion(:,:,:,2,3)); % slice, slice, 3
    Ui(:,:,:,3,ii)=imag(FFTmotion(:,:,:,2,3));
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
    
    % Check to see if there is an SPR segmentation file 
    f=[dir('*SPRregs.mat')];
    if(length(f)>0)
        copyfile(f(1).name,dname);
    end
    
    % Check for DTI file
    if(exist('DTI.mat','file'))
        copyfile('DTI.mat',dname);
    end

end

