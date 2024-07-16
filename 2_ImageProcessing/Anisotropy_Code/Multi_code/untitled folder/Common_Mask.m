function [] = Common_Mask()
set(0,'DefaultFigureWindowStyle','docked') %dock all figures
load(sprintf('/Volumes/CLJ-001/Group/drsmitty/Code/Matlab/mre_colormaps.mat'))

addpath(sprintf('%s/AP/PreInv',pwd)); 
load('t2mask.mat'); maskAP = mask; 
rmpath(sprintf('%s/AP/PreInv',pwd))

addpath(sprintf('%s/LR/PreInv',pwd)); 
load('t2mask.mat'); maskLR = mask; 
rmpath(sprintf('%s/LR/PreInv',pwd))

addpath(sprintf('%s/SI/PreInv',pwd)); 
load('t2mask.mat'); maskSI = mask; 
rmpath(sprintf('%s/SI/PreInv',pwd))

mask = maskAP.*maskLR.*maskSI;

figure;im(maskAP)
figure;im(maskLR)
figure;im(maskSI)
figure;im(mask)

cd('./AP/PostInv')
save mask.mat mask
load Mu.mat
MuNew = Mu.*mask;
save MuMasked.mat MuNew
montagestack(MuNew); caxis([0 5000]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

cd('./LR/PostInv')
save mask.mat mask
load Mu.mat
MuNew = Mu.*mask;
save MuMasked.mat MuNew
montagestack(MuNew); caxis([0 5000]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

cd('./SI/PostInv')
save mask.mat mask
load Mu.mat
MuNew = Mu.*mask;
save MuMasked.mat MuNew
montagestack(MuNew); caxis([0 5000]); set(gcf,'Colormap',stiff_color)
cd ..; cd ..

