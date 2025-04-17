% Plot hex Displacements
clear all
%close all
%nodf=getfilename('Node File Number >> ','*.nod');
%junk=load(nodf);
%xyz=junk(:,2:4);

idxf=getfilename('Index File Number >> ','*.idx');
junk=load(idxf);
idx=junk(:,2:4);

dspf=getfilename('Displacement File Number >> ','*.dsp');
junk=load(dspf);
uvw=junk(:,2:7);

dim=max(idx);
Nn=length(idx);

montsiz=input(['Size of montage (' int2str(dim(3)) ' Slices Total, Default=auto) >>']);


% Names of Displacements for figure titles
dsps(1).name='Real U Displacement';
dsps(2).name='Imag U Displacement';
dsps(3).name='Real V Displacement';
dsps(4).name='Imag V Displacement';
dsps(5).name='Real W Displacement';
dsps(6).name='Imag W Displacement';

% format tags for saved filenames
dsps(1).filetag='_ReUdsp';
dsps(2).filetag='_ImUdsp';
dsps(3).filetag='_ReVdsp';
dsps(4).filetag='_ImVdsp';
dsps(5).filetag='_ReWdsp';
dsps(6).filetag='_ImWdsp';


% Initailize dsp arrays
for ii=1:6
    dsps(ii).vals=zeros(dim);
end

for ii=1:Nn
    for jj=1:6
        dsps(jj).vals(idx(ii,1),idx(ii,2),idx(ii,3))=uvw(ii,jj);        
    end
end

for ii=1:6
    hf(ii)=figure;%(ii)
    stack=dsps(ii).vals(:,:,:);
    if(isempty(montsiz))
        montage(reshape(dsps(ii).vals(:,:,:),[dim(1) dim(2) 1 dim(3)]),'DisplayRange',[])
    else
        montage(reshape(dsps(ii).vals(:,:,:),[dim(1) dim(2) 1 dim(3)]),'DisplayRange',[],'Size',montsiz)
    end
    
    colormap(jet)
    colorbar
    title(dsps(ii).name)
    pause(0.1)
    drawnow
end
pause(0.1)
drawnow

% Save files if desired
disp('************************************************************')
disp('**                  FILE OUTPUT                           **')
disp('**          Resize images to desired size                 **')
outstm=input('**  Enter Filestem for saved files (Leave blank for no save) >>','s');

if (isempty(outstm))
    disp('****   NO IMAGES SAVED  ****')
else % Save files
    fmt=input('File Format (bmp tif jpg pdf fig eps ceps) >> ','s'); 
    for ii=1:6
        outname=[outstm dsps(ii).filetag];
        h=figure(hf(ii));
        set(gcf,'paperpositionmode','auto');
        
        if (strcmp(fmt,'ceps')) % Publication quality color eps print
            print( h, '-depsc', '-r300', outname)
        elseif (strcmp(fmt,'eps'))
            print( h, '-deps', '-r300', outname)  % Publication quality eps print
        else
            saveas(h,outname,fmt);
        end
        disp(['Figure ' int2str(ii) ' saved as ' outname '.' fmt ' in current directory'])
        clear h outname
    end
end
        

% save the results in an array
save([dspf(1:end-3) 'HexdspResults.mat'],'dsps','dspf')







