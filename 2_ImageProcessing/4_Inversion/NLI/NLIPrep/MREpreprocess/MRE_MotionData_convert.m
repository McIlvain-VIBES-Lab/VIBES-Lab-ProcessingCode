function [Ur,Ui,Motion2Image,RHcoord,voxsize_mm,freqHz,AnatomicalMRI] = MRE_MotionData_convert(MREf)
%MRE_MotionData_convert: Converts old MRE_MotionData.mat format into the
%new NLI_MRE_Input format

if(nargin<1)
    MREf='MRE_3DMotionData.mat';
end

load(MREf)
load('HeaderData.mat')

Ur=A.*cos(P);
Ui=A.*sin(P);


voxsize_mm=DirIndex(4,1:3);
rowswap=[0 1 0; 1 0 0 ; 0 0 1]; % This makes [1,:,:] = x, [:,1,:]=y, [:,:,1]=z
Motion2Image=rowswap*DirIndex(1:3,4:6);

AnatomicalMRI=MagIm;

RHcoord=DirIndex(3,3)~=-1;

save NLI_MRE_Input Ur Ui voxsize_mm RHcoord Motion2Image AnatomicalMRI freqHz

end

