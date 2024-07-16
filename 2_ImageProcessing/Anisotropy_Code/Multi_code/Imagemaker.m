%% AP image creation
cd AP
load('BMR.mat');load('SlowFastWaves.mat');load('Waveprop_BMR.mat')
XmotAP = Brfo_all(:,:,:,1); YmotAP = Brfo_all(:,:,:,2); ZmotAP = Brfo_all(:,:,:,3);
figure;imagesc(XmotAP(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
figure;imagesc(YmotAP(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
figure;imagesc(ZmotAP(:,:,30));set(gcf,'Colormap',wave_color);axis equal;
Usd_bar = flip(permute(Usd_bar,[2 1 3]),1); Ufd_bar = flip(permute(Ufd_bar,[2 1 3]),1);
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray; 
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]);axis equal; colormap hot
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray; 
figure;imagesc(abs(Ufd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray;
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_norm(:,:,30))); caxis([0 1]);axis equal; colormap hot
cd PostInv
load('Mu.mat')
figure;imagesc(Mu(:,:,30));caxis([0 5000]);set(gcf,'Colormap',stiff_color);axis equal;
cd ..
cd ..
%% LR image creation
cd LR
load('BMR.mat');load('SlowFastWaves.mat');load('Waveprop_BMR.mat')
XmotLR = Brfo_all(:,:,:,1); YmotLR = Brfo_all(:,:,:,2); ZmotLR = Brfo_all(:,:,:,3);
figure;imagesc(XmotLR(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
figure;imagesc(YmotLR(:,:,30));set(gcf,'Colormap',wave_color);axis equal;
figure;imagesc(ZmotLR(:,:,30));set(gcf,'Colormap',wave_color);axis equal;
Usd_bar = flip(permute(Usd_bar,[2 1 3]),1); Ufd_bar = flip(permute(Ufd_bar,[2 1 3]),1);
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray; 
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray;
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]);axis equal; colormap hot
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray; 
figure;imagesc(abs(Ufd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray;
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_norm(:,:,30))); caxis([0 1]);axis equal; colormap hot
cd PostInv
load('Mu.mat')
figure;imagesc(Mu(:,:,30));caxis([0 5000]);set(gcf,'Colormap',stiff_color);axis equal;
cd .. 
cd ..
%% SI image creation
cd SI
load('BMR.mat');load('SlowFastWaves.mat');load('Waveprop_BMR.mat')
XmotSI = Brfo_all(:,:,:,1); YmotSI = Brfo_all(:,:,:,2); ZmotSI = Brfo_all(:,:,:,3);
figure;imagesc(XmotSI(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
figure;imagesc(YmotSI(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
figure;imagesc(ZmotSI(:,:,30));set(gcf,'Colormap',wave_color);axis equal; 
Usd_bar = flip(permute(Usd_bar,[2 1 3]),1); Ufd_bar = flip(permute(Ufd_bar,[2 1 3]),1);
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray;
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]); axis equal; colormap gray;
figure;imagesc(abs(Usd_bar(:,:,30)));caxis([0 13]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_bar(:,:,30)));caxis([0 13]);axis equal; colormap hot;
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray; 
figure;imagesc(abs(Ufd_norm(:,:,30)));caxis([0 1]); axis equal; colormap gray;
figure;imagesc(abs(Usd_norm(:,:,30)));caxis([0 1]); axis equal; colormap hot; 
figure;imagesc(abs(Ufd_norm(:,:,30))); caxis([0 1]);axis equal; colormap hot;
cd PostInv
load('Mu.mat')
figure;imagesc(Mu(:,:,30));caxis([0 5000]);set(gcf,'Colormap',stiff_color);axis equal;
cd ..
cd ..