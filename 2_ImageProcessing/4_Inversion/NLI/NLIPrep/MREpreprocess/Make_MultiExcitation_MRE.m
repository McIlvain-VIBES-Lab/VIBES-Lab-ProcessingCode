clear all
load MRE_AP/NLI_MRE_Input.mat

Ur2(:,:,:,:,1)=Ur;
Ui2(:,:,:,:,1)=Ui;

load MRE_LR/NLI_MRE_Input.mat
Ur2(:,:,:,:,2)=Ur;
Ui2(:,:,:,:,2)=Ui;

Ur=Ur2;
Ui=Ui2;
clear Ur2 Ui2
freqHz=[freqHz freqHz];

mkdir APLR
!cp MRE_AP/* APLR
!rm APLR/NLI_MRE_Input.mat
!rm APLR/MRE_3DMotionData.mat
!rm APLR/HeaderData.mat


save('APLR/NLI_MRE_Input.mat')
