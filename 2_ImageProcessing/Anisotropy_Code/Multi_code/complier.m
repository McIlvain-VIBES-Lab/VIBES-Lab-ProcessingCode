addpath(sprintf('%s/AP/PreInv',pwd))
maskAP = load('t2mask.mat');
APmask = maskAP.mask;
rmpath(sprintf('%s/AP/PreInv',pwd))

addpath(sprintf('%s/SI/PreInv',pwd))
maskSI = load('t2mask.mat');
SImask = maskSI.mask;
rmpath(sprintf('%s/SI/PreInv',pwd))

addpath(sprintf('%s/LR/PreInv',pwd))
maskLR = load('t2mask.mat');
LRmask = maskLR.mask;
rmpath(sprintf('%s/LR/PreInv',pwd))