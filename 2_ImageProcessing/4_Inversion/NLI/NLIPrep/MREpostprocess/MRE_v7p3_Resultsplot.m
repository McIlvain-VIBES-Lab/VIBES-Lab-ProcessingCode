% Plot the distributions from a processed mat file. There is an option in MRE_plotv7p3_mask_Linear
% to supress the plots, this code just plots the processed mat files for
% easy viewing

clear svfig

[fname]=getfilename('Mat file to plot','*.mat')
load(fname)

MagIm=double(MagIm);

cmp='gray';
fsz=18;

limpercentile=99;

Imask=RealShear>0.1; % Values that are actually reconstructed

montagestack(MagIm)
colormap(gca,gray);
colorbar
title(['MR Magnitude: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(MagIm(Imask),limpercentile);
caxis([0 lim]);

H(1)=gcf;

montagestack(RealShear)
colormap(gca,cmp);
colorbar
title(['Real Shear Mod: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(RealShear(Imask),limpercentile);
caxis([0 lim]);

H(2)=gcf;

montagestack(ImagShear)
colormap(gca,cmp);
colorbar
title(['Imag Shear Mod: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(ImagShear(Imask),limpercentile);
caxis([0 lim]);

H(3)=gcf;

montagestack(DR)
colormap(gca,cmp);
colorbar
title(['Damping Ratio: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(DR(Imask),limpercentile);
caxis([0 lim]);

H(4)=gcf;

svfig=input('Format to save figures (fig, tif etc., blank for no save) >> ','s');

if(~isempty(svfig))    
    saveas(H(1),[fname(1:end-4) 'MagIm.' svfig],svfig)
    saveas(H(2),[fname(1:end-4) 'RSM.' svfig],svfig)
    saveas(H(3),[fname(1:end-4) 'ISM.' svfig],svfig)
    saveas(H(4),[fname(1:end-4) 'DR.' svfig],svfig)
end


