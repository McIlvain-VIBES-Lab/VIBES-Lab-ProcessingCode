% Plot the distributions from a processed mat file. There is an option in MRE_plotv7p3_mask_Linear
% to supress the plots, this code just plots the processed mat files for
% easy viewing

clear svfig

[fname]=getfilename('Mat file to plot','*.mat')
load(fname)

if(exist('invparam','var'))
    model=invparam.modeltype;
else
    model=1;
end

MagIm=double(MagIm);

cmp='gray';
fsz=18;

limpercentile=99;

Imask=RealShear>0.1; % Values that are actually reconstructed

montagestack(MagIm)
H(1)=gcf;
colormap(gca,gray);
colorbar
title(gca,['MR Magnitude: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(MagIm(Imask),limpercentile);
caxis(gca,[0 lim]);
pause(0.1)
drawnow




montagestack(RealShear)
H(2)=gcf;
colormap(gca,cmp);
colorbar
title(gca,['Real Shear Mod: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(RealShear(Imask),limpercentile);
caxis(gca,[0 lim]);
pause(0.1)
drawnow



montagestack(ImagShear)
H(3)=gcf;
colormap(gca,cmp);
colorbar
title(gca,['Imag Shear Mod: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(ImagShear(Imask),limpercentile);
caxis(gca,[0 lim]);
pause(0.1)
drawnow


montagestack(DR)
H(4)=gcf;
colormap(gca,cmp);
colorbar
title(gca,['Damping Ratio: ' fname],'interpreter','none')
set(gca,'fontsize',fsz)
lim=prctile(DR(Imask),limpercentile);
caxis(gca,[0 lim]);
pause(0.1)
drawnow


svfig=input('Format to save figures (fig, tif etc., blank for no save) >> ','s');

if(~isempty(svfig))    
    saveas(H(1),[fname(1:end-4) 'MagIm.' svfig],svfig)
    saveas(H(2),[fname(1:end-4) 'RSM.' svfig],svfig)
    saveas(H(3),[fname(1:end-4) 'ISM.' svfig],svfig)
    saveas(H(4),[fname(1:end-4) 'DR.' svfig],svfig)
end


