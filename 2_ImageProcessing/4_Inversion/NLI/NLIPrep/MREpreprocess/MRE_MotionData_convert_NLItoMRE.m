function [A,P,DirIndex,freqHz,MagIm] = MRE_MotionData_convert_NLItoMRE(NLIf)
% Converts new NLI_MRE_Input format to the old MRE_MotionData.mat format
%

if(nargin<1)
    NLIf='NLI_Input.mat';
end

load(NLIf)

A=abs(Ur+1i*Ui);
P=angle(Ur+1i*Ui);
DirIndex=zeros(4,6);
DirIndex(1:3,1:3)=eye(3);
DirIndex(4,1:3)=voxsize_mm;
rowswap=[0 1 0; 1 0 0 ; 0 0 1]; % This makes [1,:,:] = x, [:,1,:]=y, [:,:,1]=z
DirIndex(1:3,4:6)=rowswap'*Motion2Image;

MagIm=AnatomicalMRI;

%RHcoord=DirIndex(3,3)~=-1;

save MRE_3DMotionData A P MagIm
save HeaderData DirIndex freqHz

end