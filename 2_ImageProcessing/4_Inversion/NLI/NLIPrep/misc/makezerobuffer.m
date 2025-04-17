%% Makezerobuffer
clear all
close all

load MRE_3DMotionData
load Mask

Atmp=A;
Ptmp=P;
MagImtmp=MagIm;
masktmp=mask;
dim=size(mask);

A=zeros([dim+2 3]);
P=zeros([dim+2 3]);
mask=zeros(dim+2);
MagIm=zeros(dim+2);

A(2:dim(1)+1,2:dim(2)+1,2:dim(3)+1,:)=Atmp;
P(2:dim(1)+1,2:dim(2)+1,2:dim(3)+1,:)=Ptmp;
mask(2:dim(1)+1,2:dim(2)+1,2:dim(3)+1)=masktmp;
MagIm(2:dim(1)+1,2:dim(2)+1,2:dim(3)+1,:)=MagImtmp;

save MRE_3DMotionData A P MagIm
save Mask mask