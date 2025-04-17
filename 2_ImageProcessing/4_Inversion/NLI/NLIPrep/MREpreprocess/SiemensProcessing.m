%clear all; close all; 

files = dir; 
% Get a logical vector that tells which is a directory.
subFolders = files([files.isdir]);
subFolders(1:2)=[];
% Extract only those that are directories.
for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end

PFolder(1)=input('Select First Phase Directory: >>');
PFolder(2)=input('Select Second Phase Directory: >>');
PFolder(3)=input('Select Third Phase Directory: >>');


for ii=1:3
    cd(subFolders(PFolder(ii)).name); 
    firstfile=dir('*.IMA'); 
    firstfile=firstfile(1).name;
    DicomInfo(ii)=dicominfo(firstfile);
    if strcmp(DicomInfo(ii).ImageType(1:17),'DERIVED\PRIMARY\P')==0 && strcmp(DicomInfo(ii).ImageType(1:17),'DERIVED\PRIMARY\P\ND\RETRO')==0
        error('THESE ARE NOT PHASE IMAGES!')
    end
    if strcmp(DicomInfo(ii).Private_0051_1016(4),'P')==0 && strcmp(DicomInfo(ii).Private_0051_1016(4),'D')==0;
        error('THESE ARE NOT PHASE IMAGES (catch 2)!')
    end
    order{ii}=DicomInfo(ii).Private_0051_1014(end-1:end);
    cd ../
end

if length(unique(order))<3
    error('Three unique phase directions not found!');
end
%Checks to make sure all slice directions are the same. 
if isequal(DicomInfo(1).Private_0051_100e,DicomInfo(2).Private_0051_100e,DicomInfo(3).Private_0051_100e)==0
    error('These Phase directories do not have the same slice direction');
else 
    strcat('Slice Direction is: ',DicomInfo(1).Private_0051_100e) 
end

%Checks to make sure all pixel sizes are the same. 
if isequal(DicomInfo(1).PixelSpacing,DicomInfo(2).PixelSpacing,DicomInfo(3).PixelSpacing)==0 && ...
        isequal(DicomInfo(1).SpacingBetweenSlices,DicomInfo(2).SpacingBetweenSlices,DicomInfo(3).SpacingBetweenSlices)==0
    error('These Phase directories do not have the same voxel size');
else 
    sprintf('Voxel Size: %0.2fmm x %0.2fmm x %0.2fmm ',DicomInfo(1).PixelSpacing, DicomInfo(1).SpacingBetweenSlices) 
end

if isequal(DicomInfo(1).Rows,DicomInfo(2).Rows,DicomInfo(3).Rows)==0 && ...
        isequal(DicomInfo(1).Columns,DicomInfo(2).Columns,DicomInfo(3).Columns)==0
    error('These Phase directories do not have the same # of Rows and/or Columns');
else 
    sprintf('Matrix Size: %dx%d',DicomInfo(1).Rows, DicomInfo(1).Columns) 
end

if isequal(DicomInfo(1).RescaleSlope,DicomInfo(2).RescaleSlope,DicomInfo(3).RescaleSlope)==0 && ...
        isequal(DicomInfo(1).RescaleIntercept,DicomInfo(2).RescaleIntercept,DicomInfo(3).RescaleIntercept)==0
    error('These Phase directories do not have the same Rescale Slope/Intercept');
else 
    sprintf('Rescale Slope: %d \nRescale Intercept: %d',DicomInfo(1).RescaleSlope , DicomInfo(1).RescaleIntercept) 
end

if isequal(DicomInfo(1).CardiacNumberOfImages,DicomInfo(2).CardiacNumberOfImages,DicomInfo(3).CardiacNumberOfImages)==0
    error('These Phase directories do not have the number of cardiac cycles');
else 
    sprintf('Number of cardiac: %d',DicomInfo(1).CardiacNumberOfImages) 
end

if isequal(DicomInfo(1).Private_0051_1013,DicomInfo(2).Private_0051_1013,DicomInfo(3).Private_0051_1013)==0
    warning('The three motions directions do not have the same directions defined as "positive" motion.');
else 
    strcat('Positive Motion Direction: ', DicomInfo(1).Private_0051_1013) 
end


if strcmp(DicomInfo(1).PatientPosition,'HFS')==0 || strcmp(DicomInfo(2).PatientPosition,'HFS')==0 || strcmp(DicomInfo(3).PatientPosition,'HFS')==0
    error('You do not have a Head First Supine data set. Need constant flow motions to figure out directions');
else 
    strcat('Patient Position: ', DicomInfo(1).PatientPosition) 
end


if strcmp(DicomInfo(1).Private_0051_100e,'Cor')==1 %Coronal Slices
    index(1) = find(strcmp(order, 'rl')); %Index of RL motion
    index(2) = find(strcmp(order, 'gh')); %Index of through motion AP
    index(3) = find(strcmp(order, 'fh')); %Index of FH motion
elseif strcmp(DicomInfo(1).Private_0051_100e,'Tra')==1 %Transverse Slices
    index(1) = find(strcmp(order, 'rl')); %Index of RL motion
    index(2) = find(strcmp(order, 'ap')); %Index of AP motion
    index(3) = find(strcmp(order, 'gh')); %Index of through motion (FH)
elseif strcmp(DicomInfo(1).Private_0051_100e,'Sag')==1 %Transverse Slices
    index(1) = find(strcmp(order, 'gh')); %Index of through motion (RL) 
    index(2) = find(strcmp(order, 'ap')); %Index of AP motion
    index(3) = find(strcmp(order, 'fh')); %Index of FH
end

nslices=regexp(DicomInfo(1).Private_0051_100a,'\d*','Match'); 
nslices=str2num(nslices{1,end});
ndir=3;
venc=regexp(DicomInfo(1).Private_0051_1014,'\d*','Match'); 
venc=str2num(venc{1,end});

RawMotion=zeros(DicomInfo(1).Rows,DicomInfo(1).Columns,nslices,ndir,DicomInfo(1).CardiacNumberOfImages); 



for jj=1:3
    cd(subFolders(PFolder(index(jj))).name); %Gets RL motions first, then AP, then FH
    [dicom_images]=dir('*.IMA'); 
    lut=linspace(1,length(dicom_images),length(dicom_images));
    lut=reshape(lut,[nslices,DicomInfo(1).CardiacNumberOfImages]);
    for kk=1:length(dicom_images);
        dicom_image=dicominfo(dicom_images(kk).name); %reads in each dicom header
        [slice,phase]=find(lut==dicom_image.InstanceNumber); %finds the slice number and phase for each image based on look up table
        RawMotion(:,:,slice,jj,phase)=dicomread(dicom_images(kk).name); 
    end
    cd ../
end

FFT_Bin=input('Is this constant or harmonic motion? 1=constant, 2=Harmonic(Default)>>'); 
if isempty(FFT_Bin)==1
    FFT_Bin=2; 
elseif FFT_Bin~=1 && FFT_Bin~=2
    error('Error: Please enter 1(constant) or 2(Harmonic) motion');
end

OrientationDirection={'+LPH',
                      '+RPH',%1
                      '+LAH',%2
                      '+LPF',%3
                      '+RAH',%12
                      '+RPF',%13
                      '+LAF',%23
                      '+RAF',%123
                      '-LPH',
                      '-RPH',%1
                      '-LAH',%2
                      '-LPF',%3
                      '-RAH',%12
                      '-RPF',%13
                      '-LAF',%23
                      '-RAF',%123
                      };
OR_Index=[[1,1,1],
          [-1,1,1],
          [1,-1,1],
          [1,1,-1],
          [-1,-1,1],
          [-1,1,-1],
          [1,-1,-1],
          [-1,-1,-1],
          [-1,-1,-1],
          [1,-1,-1],
          [-1,1,-1],
          [-1,-1,1],
          [1,1,-1],
          [1,-1,1],
          [-1,1,1],
          [1,1,1]]; 
%MotionFlips=struct('OrientationDirection',OrientationDirection,'Index',OR_Index);

RawMotion=double(RawMotion);
%RawMotion=RawMotion.*DicomInfo(1).RescaleSlope+DicomInfo(1).RescaleIntercept;
RawMotion=(RawMotion-2048).*venc;
montagestack(RawMotion); colorbar

%nph=DicomInfo(1).CardiacNumberOfImages;
%nsl=size(RawMotion,3)/DicomInfo(1).CardiacNumberOfImages;
%RawMotion=reshape(RawMotion,[DicomInfo(1).Rows, DicomInfo(1).Columns,nsl,nph,3]);

%Multiplying by - if the Positive direction isn't R>>L, A>>P, F>>H (Foot to
%Head)
if strcmp(DicomInfo(1).Private_0051_1013,'+LPH')==0
    OR_IndexID=find(ismember(OrientationDirection,DicomInfo(1).Private_0051_1013)==1);
    RawMotion(:,:,:,:,1)=RawMotion(:,:,:,:,1).*OR_Index(OR_IndexID,1);
    RawMotion(:,:,:,:,2)=RawMotion(:,:,:,:,2).*OR_Index(OR_IndexID,2);
    RawMotion(:,:,:,:,3)=RawMotion(:,:,:,:,3).*OR_Index(OR_IndexID,3);
end

if strcmp(DicomInfo(1).Private_0051_1013,'+LPH')==0 || strcmp(DicomInfo(2).PatientPosition,'HFS')==0 || strcmp(DicomInfo(3).PatientPosition,'HFS')==0
    error('You do not have a Head First Supine data set. Need constant flow motions to figure out directions');
else 
    strcat('Patient Position: ', DicomInfo(1).PatientPosition) 
end



for jj=1:nslices; 
    bf_1=fft(RawMotion(:,:,jj,1,:),[],5);
    bf_2=fft(RawMotion(:,:,jj,2,:),[],5);
    bf_3=fft(RawMotion(:,:,jj,3,:),[],5);

    Motion(:,:,jj,1)=bf_1(:,:,1,1,FFT_Bin);
    Motion(:,:,jj,2)=bf_2(:,:,1,1,FFT_Bin);
    Motion(:,:,jj,3)=bf_3(:,:,1,1,FFT_Bin);
    
end
Ur=real(Motion); 
Ui=imag(Motion); 
save RawMotion.mat RawMotion Motion
%Ur=A.*cos(P); 
%Ui=A.*sin(P); 
montagestack(Ur);title('Real Displacements');colorbar; 

montagestack(Ui);title('Imag Displacements');colorbar;


voxsize_mm=[DicomInfo(1).PixelSpacing(1), DicomInfo(1).PixelSpacing(2),DicomInfo(1).SpacingBetweenSlices];
RHcoord=1; %Right handed coordinate system ensured by +LPH Motions
FOV=regexp(DicomInfo(1).Private_0051_100c,'\d*','Match');
FOV=[str2num(FOV{1}),str2num(FOV{2})];
MatSize=double([DicomInfo(1).Rows,DicomInfo(1).Columns]);
res=[FOV./MatSize,DicomInfo(1).SliceThickness];

%DOUBLE CHECK THESE!!!!!!!!!!!!
freqHz=1;
AnatomicalMRI=zeros([DicomInfo(1).Rows,DicomInfo(1).Columns,nslices]); %Can pass T2W into AnatomicalMRI variable
M2I=Motion2ImageCalc(DicomInfo(1),2);
Motion2Image=M2I;
save NLI_Input.mat  Ui Ur Motion2Image voxsize_mm freqHz AnatomicalMRI RHcoord res