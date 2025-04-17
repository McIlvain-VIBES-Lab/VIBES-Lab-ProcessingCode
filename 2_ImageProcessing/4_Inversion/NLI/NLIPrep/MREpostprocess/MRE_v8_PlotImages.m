% Plot the distributions from a processed mat file. There is an option in MRE_plotv7p3_mask_Linear
% to supress the plots, this code just plots the processed mat files for
% easy viewing

clear svfig

[fname]=getfilename('Mat file to plot','*.mat');
props=load(fname);
fields=fieldnames(props);

cmp='gray';
fsz=14;
limpercentilehigh=99;
limpercentilelow=1;

Imask=props.RealShear>0.1; % Values that are actually reconstructed

nf=0;
for ii=1:length(fields)
    test=getfield(props,fields{ii});
    if(ndims(test)==3)
        montagestack(test)
        nf=nf+1;
        H(nf)=gcf;        
        if(strcmp(fields{ii},'MagIm'))
            colormap(gca,gray);
        else
            colormap(gca,cmp);
        end
        colorbar
        drawnow
        pause(0.5)
        title(gca,[fields{ii}],'interpreter','none')        
        set(gca,'fontsize',fsz)
        pause(0.1)
        drawnow
        xlabel(['File: ' fname],'fontsize',10,'interpreter','none')
        limhigh=prctile(test(Imask),limpercentilehigh);
        limlow=prctile(test(Imask),limpercentilelow);
        caxis(gca,[min(0,limlow) max(0,limhigh)]);
        drawnow
        pause(0.5)
        
    end
end

convfigs=dir([fname(1:end-21) '*conv*.fig']);
for ii=1:length(convfigs)
    open(convfigs(ii).name);
end


svfig=input('Format to save figures (fig, tif etc., blank for no save) >> ','s');

if(~isempty(svfig))    
    for ii=1:length(H)
        saveas(H(ii),[fname(1:end-4) fields{ii} '.' svfig],svfig)
    end
end


