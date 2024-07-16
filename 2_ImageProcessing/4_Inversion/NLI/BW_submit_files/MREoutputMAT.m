function MREoutputMAT(Nplot)
% MRE_plotv7: Plots Results for MREv7 reconstructions.
% Versions:
% v1 : Original during initial development of MREv7 - very basic.
% v2 : Interpolates all values back to MR voxels for display.
% v3 : Works with elemental basis support.
% v4 : Major update for reconstuctions using MRE-Zone-v7.04+ reconstruction
%      and MRIhexmeshinterpv1 meshing step.
%       - Records reconstuction and mesh data in the same structure as the
%         results.
%       - Uses data saved from the meshing step to interpolate material
%         properties back onto the original MR voxels.





%close all;clear all;clc

%% Choose which Reconstructed Dataset to use
d=dir('*01.nod');
for ii=1:length(d)
    disp(['File ' int2str(ii) ' :: ' d(ii).name])
end
n=[];%input('Number of file to use (default = 1)  >> ');
if(isempty(n))
    n=1;
end
nodstm=d(n).name(1:end-14);

d=dir([nodstm(1:end-1) '.RE..*.prop.01.mtr']);
if length(d)~=0
    mfind=false;
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
else
    d=dir([nodstm(1:end-1) '.RE..*.prop.01.mf.mtr']);
    mfind=true;
    for ii=1:length(d)
        disp(['File ' int2str(ii) ' :: ' d(ii).name])
    end
end
clear n
n=Nplot;%input('Number of file to use (default = last)  >> ');
if(isempty(n))
    n=length(d);
end

mtrf=d(n).name;

convind=false;
avind=false;
stdind=false;
if ~mfind
    mtrstm=d(n).name(1:end-20);
    itrstr=d(n).name(end-15:end-12);
    if strcmp(mtrstm(end-4:end-1),'from')
        convind=true;
        if strcmp(mtrstm(end-6:end-5),'av')
            avind=true;
            mtrstm=d(n).name(1:end-31);
            itrstr=d(n).name(end-15:end-12);
            avstm=d(n).name(end-26:end-16);
        elseif strcmp(mtrstm(end-7:end-5),'std')
            stdind=true;
            mtrstm=d(n).name(1:end-32);
            itrstr=d(n).name(end-15:end-12);
            stdstm=d(n).name(end-27:end-16);
        end
    end
else
    mtrstm=d(n).name(1:end-23);
    itrstr=d(n).name(end-18:end-15);
    if strcmp(mtrstm(end-4:end-1),'from')
        convind=true;
        if strcmp(mtrstm(end-6:end-5),'av')
            avind=true;
            mtrstm=d(n).name(1:end-34);
            itrstr=d(n).name(end-18:end-15);
            avstm=d(n).name(end-29:end-19);
        elseif strcmp(mtrstm(end-7:end-5),'std')
            stdind=true;
            mtrstm=d(n).name(1:end-35);
            itrstr=d(n).name(end-18:end-15);
            stdstm=d(n).name(end-30:end-19);
        end
    end
end

mshindf=[nodstm 'meshind'];
basisf=[nodstm 'basis'];
meshind=load(mshindf);
rcnindf=[nodstm 'reconind'];
if(exist(rcnindf,'file'))
    fid=fopen(rcnindf);
    f=fgetl(fid);
    f=f(f~=' ');
    rcnind=f=='T';
else
    rcnind=[true true false true true false];
end

if(exist(basisf,'file'))
    basis=load(basisf);
else
    basis=[1 1 1];
end
if (sum(basis==2)>0)
    elmbasis=true;
else
    elmbasis=false;
end

filtind=[];%input('Enter <1> For Gaussian Filtering: ');
if filtind==1
    filtwind=input('Filter Window Size (default = 5): ');
    if isempty(filtwind)
        filtwind=5;
    end
    filtsig=input('Filter Strength (default = 0.05): ');
    if isempty(filtsig)
        filtsig=0.05;
    end
end

d1=dir('../*.nod');
% remove '.hom.nod' entries
ifle=0;
for ii=1:length(d1)
    if(length(d1(ii).name)>6)
        if(~strcmp(d1(ii).name(end-7:end),'.hom.nod'))
            ifle=ifle+1;
            dnod(ifle).name=d1(ii).name;
        end
    end
end

%displacement mesh
for ii=1:length(dnod)
    disp(['File ' int2str(ii) ' :: ' dnod(ii).name])
end
clear n
n=[];%input('Number of file to use (default = first)  >> ');
if(isempty(n))
%     n=length(dnod);
    n=1;
end
nod27=load(['../' dnod(n).name]);
nod27=nod27(:,2:4);
idx27=load(['../' dnod(n).name(1:end-3) 'idx']);
idx27=idx27(:,2:4);
for ii=1:3
    idx27(:,ii)=idx27(:,ii)-min(idx27(:,ii)) + 1;
end
MagImf=['../' dnod(n).name(1:end-4) 'magim.mtrin'];
if exist(MagImf,'file')
    MagIm27=load(MagImf); %Dartmouth Format
else
    MagImf=['../' dnod(n).name(1:end-4) '.magim'];
    MagIm27=load(MagImf); %UC Format
end


d=dir([nodstm 'mtrmesh*nod']);
nmesh=length(d);
for ii=1:nmesh
    junk=load(d(ii).name);
    mesh(ii).nod(:,:)=junk(:,2:4);
    
    if(exist([d(ii).name(1:end-3) 'elm'],'file'))
        junk=load([d(ii).name(1:end-3) 'elm']);
        mesh(ii).in=junk(:,2:9);
        mesh(ii).ne=length(junk);
    elseif(elmbasis)
        error('Elemental basis requires element file')
    end
    
    if(elmbasis)
        % Find element centroids to locate material property points
        for jj=1:mesh(ii).ne
            for kk=1:3
                mesh(ii).elmcent(jj,kk)=mean( mesh(ii).nod(mesh(ii).in(jj,:),kk) );
            end
        end
    end
    
    
    % Figure out mesh resolution
    for jj=1:3
        
        unq=unique(mesh(ii).nod(:,jj));
        srt=sort(unq);
        mesh(ii).res(jj)=(srt(end)-srt(1))/(length(srt)-1);
        
        if ii==1 %consider displacement mesh
            unqd=unique(nod27(:,jj));
            srtd=sort(unqd);
        end
        
        %disp('XXX NEED TO CONSIDER BASIS HERE')
        if(jj==1)
            mesh(ii).nodX=srt;
            if ii==1
                nodX27=srtd;
            end
        elseif(jj==2)
            mesh(ii).nodY=srt;
            if ii==1
                nodY27=srtd;
            end
        elseif(jj==3)
            mesh(ii).nodZ=srt;
            if ii==1
                nodZ27=srtd;
            end
        end
        
        if(elmbasis)
            unq=unique(mesh(ii).elmcent(:,jj));
            srt=sort(unq);
            if(jj==1)
                mesh(ii).elmX=srt;
            elseif(jj==2)
                mesh(ii).elmY=srt;
            elseif(jj==3)
                mesh(ii).elmZ=srt;
            end
        end
        
        mesh(ii).nod(:,jj)=mesh(ii).nod(:,jj)-min(mesh(ii).nod(:,jj));
        
        if(elmbasis)
            mesh(ii).elmcent(:,jj)=mesh(ii).elmcent(:,jj)-min(mesh(ii).elmcent(:,jj));
        end
        
        mesh(ii).ind(:,jj)=round(mesh(ii).nod(:,jj)/mesh(ii).res(jj));
        mesh(ii).ind(:,jj)=mesh(ii).ind(:,jj) - min(mesh(ii).ind(:,jj)) + 1;
        
        if(elmbasis)
            mesh(ii).indelm(:,jj)=round((mesh(ii).elmcent(:,jj)- min(mesh(ii).elmcent(:,jj)))/mesh(ii).res(jj));
            mesh(ii).indelm(:,jj)=mesh(ii).indelm(:,jj)+ 1;
        end
    end
end






%% Interpolate all properties back onto the MR voxels for display
%VI = INTERP3(X,Y,Z,V,XI,YI,ZI)
% X,Y,Z values for original mesh from nod27.

% Define whether to interpolate to the original MR voxels (1) or the
% interpolated mesh locations (2)
%interptyp=1;
interptyp=1;%input('Interploation Type [1 = MR Voxels, 2 = Mesh Locations] (default = 1): ');
if isempty(interptyp)
    interptyp=1;
end

% THese are the mesh resolution, which may be interpolated from the
% original data
if(interptyp==1)
    disp('Interpolating Properties to Original MRE data resolution..')
    if(exist(['../' dnod(n).name(1:end-3) 'InterpLocations.mat'],'file'))
        load(['../' dnod(n).name(1:end-3) 'InterpLocations.mat'])
    else
        disp(['Interpolation File Not Found:',['../' dnod(n).name(1:end-3) 'InterpLocations.mat']])
        error(['Requires ' dnod(n).name(1:end-3) 'InterpLocations.mat' ' file: Use MREhexmesh_interpv2 or higher to generate mesh'])
    end
    Xmr=xin;
    Ymr=yin;
    Zmr=zin;
    
    %Get Mask File
    dm=dir('../*.mask.mat');
    for ii=1:length(dm)
        disp(['File ' int2str(ii) ' :: ' dm(ii).name])
    end
    n=1;%input('Number of file to use (default = 1)  >> ');
    if(isempty(n))
        n=1;
    end
    mskf=['../' dm(n).name];
    load(mskf);
%     disp(['Loading Mask and MRE3DMotionData.mat File:  ' '../' msk])
%     mskf=['../' msk];
%     if exist(mskf,'file')
%         load(mskf);
%     else
%         mskf=input('Default Mask File Not Found, Please Enter Path and Name:  >>','s');
%         load(mskf);
%     end
    % Get MagIm Values from MRE Data
    if exist(['../../../MRE_3DMotionData.mat'],'file')
        meshmag=false;
        load(['../../../MRE_3DMotionData.mat'])
    else
        % Get MagIm Values from mesh files
        meshmag=true;
        MagIm=zeros(max(idx27));
        for ii=1:length(nod27)
            MagIm(idx27(ii,1),idx27(ii,2),idx27(ii,3))=MagIm27(ii,2);
        end
    end
    
    
    
elseif(interptyp==2)
    disp('Interpolating Properties to MRE mesh resolution..')
    Xmr=sort(unique(nod27(:,1)));
    Ymr=sort(unique(nod27(:,2)));
    Zmr=sort(unique(nod27(:,3)));
    
    % Get MagIm Values from mesh files
    meshmag=true;
    MagIm=zeros(max(idx27));
    for ii=1:length(nod27)
        MagIm(idx27(ii,1),idx27(ii,2),idx27(ii,3))=MagIm27(ii,2);
    end
    mask=MagIm~=0;
end


outputdef=1;
outputind=1;%input('Enter <1> to output solution variables (default = 1): ');
if isempty(outputind)
    outputind=outputdef;
end
%% Plot properties

displayinterp=1; % Factor to interpolate images for display (set to 2 to double the size, sometimes the images are too small on screen to see.)

sliceind=[];%input('Enter Slice # for single slice plots: ');
if(isempty(sliceind))
    sliceind=0;
end
if sliceind>size(MagIm,3)
    sliceind=size(MagIm,3);
end


if ~meshmag
    clear intstack sol
    sol=MagIm(:,:,:);
    for jj=1:size(MagIm,3)
        if(displayinterp==1)
            intstack(:,:,jj)=MagIm(:,:,jj);
        else
            intstack(:,:,jj)=interp2(MagIm(:,:,jj),displayinterp);
        end
    end
    save('intstack.mat','intstack');
    if outputind==1
        save([mtrf,'.magim.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title('MR Magnitude image')
else
    MagImStack=interp3(nodY27,nodX27',nodZ27,MagIm,Ymr,Xmr',Zmr,'linear',0).*mask;
    clear intstack sol
    sol=MagImStack(:,:,:);
    for jj=1:size(MagImStack,3)
        if(displayinterp==1)
            intstack(:,:,jj)=MagImStack(:,:,jj);
        else
            intstack(:,:,jj)=interp2(MagImStack(:,:,jj),displayinterp);
        end
    end
    save('intstack.mat','intstack');
    if outputind==1
        save([mtrf,'.magim.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title('MR Magnitude image')
end

if ~convind
    d=dir([mtrstm 'RE..' itrstr '.prop.*.mtr']);
elseif avind
    d=dir([mtrstm 'RE..' avstm itrstr '.prop.*.mtr']);
elseif stdind
    d=dir([mtrstm 'RE..' stdstm itrstr '.prop.*.mtr']);
end
numprop=length(d);
mfpropind=false(numprop,1);
for ii=1:numprop
    rmsh=meshind(1,ii);
    imsh=meshind(2,ii);
    if mfind
        if d(ii).name(end-5:end-4)=='mf'
            mfpropind(ii)=true;
        end
    end
    if ~mfpropind(ii)
        if ~convind
            prop(ii).revals=load([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            mtrf=[mtrstm '.' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        elseif avind
            prop(ii).revals=load([mtrstm 'RE..' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            mtrf=[mtrstm '.' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        elseif stdind
            prop(ii).revals=load([mtrstm 'RE..' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
            mtrf=[mtrstm '.' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        end
    else
        if ~convind
            prop(ii).revals=load([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            mtrf=[mtrstm '.' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        elseif avind
            prop(ii).revals=load([mtrstm 'RE..' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            mtrf=[mtrstm '.' avstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        elseif stdind
            prop(ii).revals=load([mtrstm 'RE..' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            prop(ii).imvals=load([mtrstm 'IM..' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
            mtrf=[mtrstm '.' stdstm itrstr '.prop.' sprintf('%2.2i',ii) '.mtr'];
        end
    end
    prop(ii).npr=length(prop(ii).revals);
    prop(ii).npi=length(prop(ii).imvals);
    
    if(basis(ii)==1)
        prop(ii).restack=zeros(max(mesh(rmsh).ind));
        prop(ii).imstack=zeros(max(mesh(imsh).ind));
    elseif(basis(ii)==2)
        prop(ii).restack=zeros(max(mesh(rmsh).indelm));
        prop(ii).imstack=zeros(max(mesh(imsh).indelm));
    end
    %disp('XX Up to here, fix this.')
    for jj=1:prop(ii).npr
        if(basis(ii)==1)
            i1=mesh(rmsh).ind(jj,1);
            i2=mesh(rmsh).ind(jj,2);
            i3=mesh(rmsh).ind(jj,3);
        elseif(basis(ii)==2)
            i1=mesh(rmsh).indelm(jj,1);
            i2=mesh(rmsh).indelm(jj,2);
            i3=mesh(rmsh).indelm(jj,3);
        end
        prop(ii).restack(i1,i2,i3)=prop(ii).revals(jj,2);
        if mfpropind(ii)
            prop(ii).restack(i1,i2,i3,2)=prop(ii).revals(jj,3);
        end
    end
    for jj=1:prop(ii).npi
        if(basis(ii)==1)
            i1=mesh(imsh).ind(jj,1);
            i2=mesh(imsh).ind(jj,2);
            i3=mesh(imsh).ind(jj,3);
        elseif(basis(ii)==2)
            i1=mesh(imsh).indelm(jj,1);
            i2=mesh(imsh).indelm(jj,2);
            i3=mesh(imsh).indelm(jj,3);
        end
        prop(ii).imstack(i1,i2,i3)=prop(ii).imvals(jj,2);
        if mfpropind(ii)
            prop(ii).imstack(i1,i2,i3,2)=prop(ii).imvals(jj,3);
        end
    end
    
    % Interpolate to Chosen Resolution (MR resolution, interptyp=1, or mesh resolution, interptyp=2)
    if(basis(ii)==1)
        prop(ii).reMRstack=interp3(mesh(rmsh).nodY,mesh(rmsh).nodX',mesh(rmsh).nodZ,prop(ii).restack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        prop(ii).imMRstack=interp3(mesh(imsh).nodY,mesh(imsh).nodX',mesh(imsh).nodZ,prop(ii).imstack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        if mfpropind(ii)
            prop(ii).reMRstack(:,:,:,2)=interp3(mesh(rmsh).nodY,mesh(rmsh).nodX',mesh(rmsh).nodZ,prop(ii).restack(:,:,:,2),Ymr,Xmr',Zmr,'linear',0).*mask;
            prop(ii).imMRstack(:,:,:,2)=interp3(mesh(imsh).nodY,mesh(imsh).nodX',mesh(imsh).nodZ,prop(ii).imstack(:,:,:,2),Ymr,Xmr',Zmr,'linear',0).*mask;
        end
    elseif(basis(ii)==2)
        prop(ii).reMRstack=interp3(mesh(rmsh).elmY,mesh(rmsh).elmX',mesh(rmsh).elmZ,prop(ii).restack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        prop(ii).imMRstack=interp3(mesh(imsh).elmY,mesh(imsh).elmX',mesh(imsh).elmZ,prop(ii).imstack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        if mfpropind(ii)
            prop(ii).reMRstack(:,:,:,2)=interp3(mesh(rmsh).elmY,mesh(rmsh).elmX',mesh(rmsh).elmZ,prop(ii).restack(:,:,:,2),Ymr,Xmr',Zmr,'linear',0).*mask;
            prop(ii).imMRstack(:,:,:,2)=interp3(mesh(imsh).elmY,mesh(imsh).elmX',mesh(imsh).elmZ,prop(ii).imstack(:,:,:,2),Ymr,Xmr',Zmr,'linear',0).*mask;
        end
    end
    
    if(rcnind(2*(ii-1)+1))
        clear intstack sol
        sol=prop(ii).reMRstack(:,:,:,1);
        for jj=1:size(prop(ii).reMRstack(:,:,:,1),3)
            if(displayinterp==1)
                intstack(:,:,jj)=prop(ii).reMRstack(:,:,jj,1);
            else
                intstack(:,:,jj)=interp2(prop(ii).reMRstack(:,:,jj,1),displayinterp);
            end
        end
        if filtind==1
            intstack=stackfilt(intstack,filtwind,filtsig);
        end
        if outputind==1
            save([mtrf,'.real.mat'],'intstack');
        end
        if sliceind==0
            montagestack(intstack);
        else
            figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
        end
        colorbar
        title(['Real Prop ' int2str(ii) '; res = ' num2str(mesh(rmsh).res)])
        if mfpropind(ii)
            clear intstack sol
            sol=prop(ii).reMRstack(:,:,:,2);
            for jj=1:size(prop(ii).reMRstack(:,:,:,2),3)
                if(displayinterp==1)
                    intstack(:,:,jj)=prop(ii).reMRstack(:,:,jj,2);
                else
                    intstack(:,:,jj)=interp2(prop(ii).reMRstack(:,:,jj,2),displayinterp);
                end
            end
            if filtind==1
                intstack=stackfilt(intstack,filtwind,filtsig);
            end
            if outputind==1
                save([mtrf,'.real.mf.mat'],'intstack');
            end
            if sliceind==0
                montagestack(intstack);
            else
                figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
            end
            colorbar
            title(['Real Prop ' int2str(ii) ' \alpha component; res = ' num2str(mesh(rmsh).res)])
        end
    end
    
    if(rcnind(2*(ii-1)+2))
        clear intstack sol
        sol=prop(ii).imMRstack(:,:,:,1);
        for jj=1:size(prop(ii).imMRstack(:,:,:,1),3)
            if(displayinterp==1)
                intstack(:,:,jj)=prop(ii).imMRstack(:,:,jj,1);
            else
                intstack(:,:,jj)=interp2(prop(ii).imMRstack(:,:,jj,1),displayinterp);
            end
        end
        if filtind==1
            intstack=stackfilt(intstack,filtwind,filtsig);
        end
        if outputind==1
            save([mtrf,'.imag.mat'],'intstack');
        end
        if sliceind==0
            montagestack(intstack);
        else
            figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
        end
        colorbar
        title(['Imag Prop ' int2str(ii) '; res = ' num2str(mesh(imsh).res)])
        if mfpropind(ii)
            clear intstack sol
            sol=prop(ii).imMRstack(:,:,:,2);
            for jj=1:size(prop(ii).imMRstack(:,:,:,2),3)
                if(displayinterp==1)
                    intstack(:,:,jj)=prop(ii).imMRstack(:,:,jj,2);
                else
                    intstack(:,:,jj)=interp2(prop(ii).imMRstack(:,:,jj,2),displayinterp);
                end
            end
            if filtind==1
                intstack=stackfilt(intstack,filtwind,filtsig);
            end
            if outputind==1
                save([mtrf,'.imag.mf.mat'],'intstack');
            end
            if sliceind==0
                montagestack(intstack);
            else
                figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
            end
            colorbar
            title(['Imag Prop ' int2str(ii) ' \alpha component; res = ' num2str(mesh(imsh).res)])
        end
    end
    
    
    
    %     clear intstack
    %     for jj=1:size(prop(ii).restack,3)
    %       intstack(:,:,jj)=interp2(prop(ii).restack(:,:,jj),3);
    %     end
    %     montagestack(intstack)
    %     colorbar
    %     title(['Real Prop ' int2str(ii) ' res = ' num2str(mesh(rmsh).res)])
    %     clear intstack
    %     for jj=1:size(prop(ii).imstack,3)
    %       intstack(:,:,jj)=interp2(prop(ii).imstack(:,:,jj),2);
    %     end
    %     if(ii~=3)
    %     montagestack(intstack)
    %     colorbar
    %     title(['Imag Prop ' int2str(ii) ' res = ' num2str(mesh(imsh).res)])
    %     end
    %
    
end
%%


%%

% Build damping ratio images. Need to interpolate each property onto a
% common mesh/grid.

dfreq=50;%input('Frequency for Damping Property Display: ');

if ~mfpropind(1)
    murat=prop(1).imMRstack(:,:,:,1)./prop(1).reMRstack(:,:,:,1);
else
    murat=(prop(1).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).imMRstack(:,:,:,2)))./(prop(1).reMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).reMRstack(:,:,:,2)));
    mure=(prop(1).reMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).reMRstack(:,:,:,2)));
    muim=(prop(1).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).imMRstack(:,:,:,2)));
    dispmure=100*(prop(1).reMRstack(:,:,:,2).*(prop(1).reMRstack(:,:,:,1).*((2*pi*dfreq).^(prop(1).reMRstack(:,:,:,2)-ones(size(prop(1).reMRstack(:,:,:,2)))))));
    dispmuim=100*(prop(1).imMRstack(:,:,:,2).*(prop(1).imMRstack(:,:,:,1).*((2*pi*dfreq).^(prop(1).imMRstack(:,:,:,2)-ones(size(prop(1).imMRstack(:,:,:,2)))))));
end
if ~mfpropind(2)
    if rcnind(4)
        rhorat=prop(2).imMRstack(:,:,:,1)./prop(2).reMRstack(:,:,:,1);
    else
        rhorat=0;
    end
else
    freqinv=input('Enter <1> for Inverse Frequency Relationship for \rho_I: ');
    if isempty(freqinv)
        freqinv=0;
    end
    if freqinv==0
        rhorat=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)))./prop(2).reMRstack(:,:,:,1);
        rhoim=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)));
        disprhoim=100*(prop(2).imMRstack(:,:,:,2).*(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^(prop(2).imMRstack(:,:,:,2)-ones(size(prop(2).imMRstack(:,:,:,2)))))));
    elseif freqinv==1
        rhorat=(prop(2).imMRstack(:,:,:,1).*((1/(2*pi*dfreq)).^prop(2).imMRstack(:,:,:,2)))./prop(2).reMRstack(:,:,:,1);
        rhoim=(prop(2).imMRstack(:,:,:,1).*((1/(2*pi*dfreq)).^prop(2).imMRstack(:,:,:,2)));
        disprhoim=-100*(prop(2).imMRstack(:,:,:,2).*(prop(2).imMRstack(:,:,:,1).*((1/(2*pi*dfreq)).^(prop(2).imMRstack(:,:,:,2)+ones(size(prop(2).imMRstack(:,:,:,2)))))));
    end
end
DR=0.5.*(murat-rhorat);
DR(isnan(DR))=0.0;
DR(isinf(DR))=0.0;

if mfpropind(1)&&mfpropind(2)
    dispDR=100*((1./mure).*dispmuim-murat.*(1./mure).*dispmure-(1./rhoim).*disprhoim);
end

RC=0.5*murat./DR;
RC(isnan(RC))=0.0;
RC(isinf(RC))=0.0;

if mfpropind(1)
    clear intstack sol
    sol=mure(:,:,:);
    for jj=1:size(mure,3)
        if(displayinterp==1)
            intstack(:,:,jj)=mure(:,:,jj);
        else
            intstack(:,:,jj)=interp2(mure(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.mu_real.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\mu_R at ',num2str(dfreq),' Hz']);
    colorbar
    
    clear intstack sol
    sol=dispmure(:,:,:);
    for jj=1:size(dispmure,3)
        if(displayinterp==1)
            intstack(:,:,jj)=dispmure(:,:,jj);
        else
            intstack(:,:,jj)=interp2(dispmure(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.disp.mu_real.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\mu_R Dispersion at ',num2str(dfreq),' Hz']);
    colorbar
    
    clear intstack sol
    sol=muim(:,:,:);
    for jj=1:size(muim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=muim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(muim(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.mu_imag.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\mu_I at ',num2str(dfreq),' Hz']);
    colorbar
    
    clear intstack sol
    sol=dispmuim(:,:,:);
    for jj=1:size(dispmuim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=dispmuim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(dispmuim(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.disp.mu_imag.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\mu_I Dispersion at ',num2str(dfreq),' Hz']);
    colorbar
end

if mfpropind(2)
    clear intstack sol
    sol=rhoim(:,:,:);
    for jj=1:size(rhoim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=rhoim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(rhoim(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.rho_imag.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\rho_I at ',num2str(dfreq),' Hz']);
    colorbar
    
    clear intstack sol
    sol=disprhoim(:,:,:);
    for jj=1:size(disprhoim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=disprhoim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(disprhoim(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.disp.rho_imag.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['\rho_I Dispersion at ',num2str(dfreq),' Hz']);
    colorbar
end

clear intstack sol
sol=DR(:,:,:);
for jj=1:size(DR,3)
    if(displayinterp==1)
        intstack(:,:,jj)=DR(:,:,jj);
    else
        intstack(:,:,jj)=interp2(DR(:,:,jj),displayinterp);
    end
end
if filtind==1
    intstack=stackfilt(intstack,filtwind,filtsig);
end
if outputind==1
    save([mtrf,'.DR.at.',num2str(dfreq),'.mat'],'intstack');
end
if sliceind==0
    montagestack(intstack);
else
    figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
end
if ~mfpropind(1)
    title('Damping Ratio')
else
    title(['Damping Ratio at ',num2str(dfreq),' Hz']);
end
colorbar

if mfpropind(1)&&mfpropind(2)
    clear intstack sol
    sol=dispDR(:,:,:);
    for jj=1:size(dispDR,3)
        if(displayinterp==1)
            intstack(:,:,jj)=dispDR(:,:,jj);
        else
            intstack(:,:,jj)=interp2(dispDR(:,:,jj),displayinterp);
        end
    end
    if filtind==1
        intstack=stackfilt(intstack,filtwind,filtsig);
    end
    if outputind==1
        save([mtrf,'.disp.DR.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['Damping Ratio Dispersion at ',num2str(dfreq),' Hz']);
    colorbar
end

clear intstack sol
sol=RC(:,:,:);
for jj=1:size(RC,3)
    if(displayinterp==1)
        intstack(:,:,jj)=RC(:,:,jj);
    else
        intstack(:,:,jj)=interp2(RC(:,:,jj),displayinterp);
    end
end
if filtind==1
    intstack=stackfilt(intstack,filtwind,filtsig);
end
if outputind==1
    save([mtrf,'.RC.at.',num2str(dfreq),'.mat'],'intstack');
end
if sliceind==0
    montagestack(intstack);
else
    figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
end
if ~mfpropind(1)
    title('Rayleigh Composition')
else
    title(['Rayleigh Composition at ',num2str(dfreq),' Hz']);
end
colorbar

KKanalysisdef=0;
KKanalysis=[];%input('Enter <1> for Kramers-Kronig Analysis (default 0): ');
if isempty(KKanalysis)
    KKanalysis=KKanalysisdef;
end
if KKanalysis==1
    freqrange=input('Frequency Range (Hz) for Kramers-Kronig Analysis: ');
    
    Tomeg=20;
    omegi=2*pi*freqrange(1);
    omegf=2*pi*freqrange(2);
    domeg=(omegf-omegi)/Tomeg;
    omeg0=omegi+0.5*domeg;
    
    if ~mfpropind(1)
        mure=prop(1).reMRstack(:,:,:,1);
        muim=prop(1).imMRstack(:,:,:,1);
        rhore=prop(2).reMRstack(:,:,:,1);
        rhoim=prop(2).imMRstack(:,:,:,1);
    else
        mure=(prop(1).reMRstack(:,:,:,1).*((omegf).^prop(1).reMRstack(:,:,:,2)));
        muim=(prop(1).imMRstack(:,:,:,1).*((omegf).^prop(1).imMRstack(:,:,:,2)));
        rhore=prop(2).reMRstack(:,:,:,1);
        if freqinv~=1
            rhoim=(prop(2).imMRstack(:,:,:,1).*((omegf).^prop(2).imMRstack(:,:,:,2)));
        else
            rhoim=(prop(2).imMRstack(:,:,:,1).*((1/omegf).^prop(2).imMRstack(:,:,:,2)));
        end
    end
    beta=rhoim./rhore;
    
    K1true=(mure+beta.*muim)./(mure.^2+muim.^2);
    
    if ~mfpropind(1)
        mure=prop(1).reMRstack(:,:,:,1);
        muim=prop(1).imMRstack(:,:,:,1);
        rhore=prop(2).reMRstack(:,:,:,1);
        rhoim=prop(2).imMRstack(:,:,:,1);
    else
        mure=(prop(1).reMRstack(:,:,:,1).*((omeg0).^prop(1).reMRstack(:,:,:,2)));
        muim=(prop(1).imMRstack(:,:,:,1).*((omeg0).^prop(1).imMRstack(:,:,:,2)));
        rhore=prop(2).reMRstack(:,:,:,1);
        if freqinv~=1
            rhoim=(prop(2).imMRstack(:,:,:,1).*((omeg0).^prop(2).imMRstack(:,:,:,2)));
        else
            rhoim=(prop(2).imMRstack(:,:,:,1).*((1/omeg0).^prop(2).imMRstack(:,:,:,2)));
        end
    end
    beta=rhoim./rhore;
    
    K10=(domeg/omeg0)*(mure+beta.*muim)./(mure.^2+muim.^2);
    omegi=omeg0;
    for ii=1:Tomeg-1
        omegi=omegi+domeg;
        if ~mfpropind(1)
            mure=prop(1).reMRstack(:,:,:,1);
            muim=prop(1).imMRstack(:,:,:,1);
            rhore=prop(2).reMRstack(:,:,:,1);
            rhoim=prop(2).imMRstack(:,:,:,1);
        else
            mure=(prop(1).reMRstack(:,:,:,1).*((omegi).^prop(1).reMRstack(:,:,:,2)));
            muim=(prop(1).imMRstack(:,:,:,1).*((omegi).^prop(1).imMRstack(:,:,:,2)));
            rhore=prop(2).reMRstack(:,:,:,1);
            rhoim=(prop(2).imMRstack(:,:,:,1).*((omegi).^prop(2).imMRstack(:,:,:,2)));
            if freqinv~=1
                rhoim=(prop(2).imMRstack(:,:,:,1).*((omegi).^prop(2).imMRstack(:,:,:,2)));
            else
                rhoim=(prop(2).imMRstack(:,:,:,1).*((1/omegi).^prop(2).imMRstack(:,:,:,2)));
            end
        end
        beta=rhoim./rhore;
        K10=K10-(domeg/omegi)*(2/pi)*(beta.*mure-muim)./(mure.^2+muim.^2);
    end
    
    clear intstack sol
    sol=K1true(:,:,:);
    for jj=1:size(mure,3)
        if(displayinterp==1)
            intstack(:,:,jj)=K1true(:,:,jj);
        else
            intstack(:,:,jj)=interp2(K1true(:,:,jj),displayinterp);
        end
    end
    if outputind==1
        save([mtrf,'.K1true.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['True K_1 at Freq ',num2str(freqrange(2)),' Hz']);
    colorbar
    
    clear intstack sol
    sol=K10(:,:,:);
    for jj=1:size(K10,3)
        if(displayinterp==1)
            intstack(:,:,jj)=K10(:,:,jj);
        else
            intstack(:,:,jj)=interp2(K10(:,:,jj),displayinterp);
        end
    end
    if outputind==1
        save([mtrf,'.K1kk.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['Kramers-Kronig K_1 at Freq ',num2str(freqrange(2)),' Hz']);
    colorbar
    
    diffKK=100.0*(K10-K1true)./K1true;
    clear intstack sol
    sol=diffKK(:,:,:);
    for jj=1:size(diffKK,3)
        if(displayinterp==1)
            intstack(:,:,jj)=diffKK(:,:,jj);
        else
            intstack(:,:,jj)=interp2(diffKK(:,:,jj),displayinterp);
        end
    end
    if outputind==1
        save([mtrf,'.K1diff.at.',num2str(dfreq),'.mat'],'intstack');
    end
    if sliceind==0
        montagestack(intstack);
    else
        figure;  imagesc(intstack(:,:,sliceind)); colormap('gray');
    end
    title(['Percentage Difference Kramers-Kronig vs. True K_1 at Freq ',num2str(freqrange(2)),' Hz']);
    colorbar
    
end


%% Read and process .invparam file
invparamf=[nodstm 'invparam'];
if(exist(invparamf,'file'))
    disp('Inversion Paraters stored in invparam structure')
    fid=fopen(invparamf);
    endoffile=false;
    while(~endoffile)
        junkheader=fgetl(fid);
        vname=fgetl(fid);
        values=fgetl(fid);
        if(ischar(values)) % end of file not reached
            values=strtrim(values);
            if((length(values)>1)&&strcmp(values(1:2),'-s')) % Character value
                eval(['invparam.' strtrim(vname) ' = ''' values(3:end) ''';'])
            else % Numeric value
                eval(['invparam.' strtrim(vname) ' = [' values '];'])
            end
        else
            endoffile=true;
        end
        
    end
    if(exist(['../' invparam.nodefile(1:end-7) 'meshinput.mat'],'file'))
        meshinputs=load(['../' invparam.nodefile(1:end-7) 'meshinput.mat']);
    else
        %Get Mesh File
        dm=dir('../*.meshinput.mat');
        for ii=1:length(dm)
            disp(['File ' int2str(ii) ' :: ' dm(ii).name])
        end
        n=[];%input('Number of file to use (default = 1)  >> ');
        if(isempty(n))
            n=1;
        end
        meshf=['../' dm(n).name];
        if exist(meshf,'file')
            meshinputs=load(meshf);
        else
            error(['Requires ' '../' invparam.nodefile(1:end-7) 'meshinput.mat' ' file: Use MREhexmesh_interpv2 or higher to generate mesh'])
        end
    end
    
    fclose(fid);
else
    invparam='Data Not Present - MREv7.04 and above required for recording of inversion parameters';
end



if(exist(invparamf,'file')) % Info about the material model is present
    modeltyp=invparam.modeltype;
else % Assume it is model = 1 (iso incompressible, all versions where invparamf is not present can only handle isocompressible
    modeltyp=1;
end

%% Save properties
RealShear=prop(1).reMRstack;
ImagShear=prop(1).imMRstack;
RealDens=prop(2).reMRstack;
ImagDens=prop(2).imMRstack;
if(modeltyp==1) % iso incompressible
    RealBulk=prop(3).reMRstack;
    ImagBulk=prop(3).imMRstack;
elseif(modeltyp==3) % iso incompressible
    RealLambda=prop(3).reMRstack;
    ImagLambda=prop(3).imMRstack;
end


for ii=1:nmesh
    meshres(ii,:)=mesh(ii).res;
end



outmatf=[mtrstm '.' itrstr '.ReconProps.mat'];
recondir=pwd;
save(outmatf,'MagIm','RealShear','ImagShear','RealDens','ImagDens','DR','RC','meshind','meshres','recondir');
if(exist(invparamf,'file')) % Info about the inversion params is present
    save(outmatf,'meshinputs','invparam','-append')
end
if(modeltyp==1) % iso incompressible
    save(outmatf,'RealBulk','ImagBulk','-append')
elseif(modeltyp==3) % iso incompressible
    save(outmatf,'RealLambda','ImagLambda','-append')
end

















