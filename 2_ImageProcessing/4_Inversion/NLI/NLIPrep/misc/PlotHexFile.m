% Plot hex Displacements
clear all
close all
%nodf=getfilename('Node File Number >> ','*.nod');
%junk=load(nodf);
%xyz=junk(:,2:4);

nodf=getfilename('Node File Number >> ','*.nod');
junk=load(nodf);
nod=junk(:,2:4);

% ROund nod to 8 s.f. to avoid rounding errors
nod=round(nod,8,'significant');

hexf=getfilename('Output File Number >> ','*');
junk=load(hexf);
nval=size(junk,2)-1;
uvw=junk(:,2:end);

disp(['File : ' hexf])

idx=zeros(size(nod));
res=zeros(1,3);
% Compute index Locations from node file 
for ii=1:3
    nod(:,ii)=nod(:,ii)-min(nod(:,ii));
    res(ii)=mean(diff(unique(sort(nod(:,ii)))));
    idx(:,ii)=round(nod(:,ii)./res(ii))+1;
    %idx(:,ii)=idx(:,ii)-min(idx(:,ii))+1;    
end
dim=max(idx);
Nn=length(idx);

montsiz=input(['Size of montage (' int2str(dim(3)) ' Slices Total, Default=auto) >>']);

% Initailize dsp arrays
for ii=1:nval
    dsps(ii).name=['Hex val column ' int2str(ii+1)];
    dsps(ii).filetag=['Hex_col' int2str(ii+1)];
    dsps(ii).vals=zeros(dim);
end

for ii=1:Nn
    for jj=1:nval
        dsps(jj).vals(idx(ii,1),idx(ii,2),idx(ii,3))=uvw(ii,jj);        
    end
end

for ii=1:nval
    figure(ii)
    stack=dsps(ii).vals(:,:,:);
    if(isempty(montsiz))
        montage(reshape(dsps(ii).vals(:,:,:),[dim(1) dim(2) 1 dim(3)]),'DisplayRange',[])
    else
        montage(reshape(dsps(ii).vals(:,:,:),[dim(1) dim(2) 1 dim(3)]),'DisplayRange',[],'Size',montsiz)
    end
    
    colormap(jet)
    colorbar
    title(dsps(ii).name)
    drawnow
    pause(0.1)
end

% Save files if desired
disp('************************************************************')
disp('**                  FILE OUTPUT                           **')
disp('**          Resize images to desired size                 **')
outstm=input('**  Enter Filestem for saved files (Leave blank for no save) >>','s');

if (isempty(outstm))
    disp('****   NO IMAGES SAVED  ****')
else % Save files
    fmt=input('File Format (bmp tif jpg pdf fig eps ceps) >> ','s'); 
    for ii=1:nval
        outname=[outstm dsps(ii).filetag];
        h=figure(ii);
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
vals=dsps;
save('HexFileResults.mat','vals','hexf')







