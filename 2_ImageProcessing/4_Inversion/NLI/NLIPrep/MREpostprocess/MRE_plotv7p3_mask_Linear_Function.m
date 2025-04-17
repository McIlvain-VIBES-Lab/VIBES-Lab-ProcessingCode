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

function MRE_plotv7p3_mask_Linear_Function(nodf,mtrf,showimg)


if(nargin<3) % Default is to not show images 
  showimg=false;
end

if(nargin<2)
    error('must supply material nod file and mtr filename')
end


%% Reconstructed Dataset to use

nodstm=nodf(1:end-14);

if(strcmp(mtrf(end-6:end),'mf.mtr'))
    mfind=true;
else
    mfind=false;
end

if ~mfind
    mtrstm=mtrf(1:end-20);
    itrstr=mtrf(end-15:end-12);
else
    mtrstm=mtrf(1:end-23);
    itrstr=mtrf(end-18:end-15);
end    
    


% Define whether to interpolate to the original MR voxels (1) or the
% interpolated mesh locations (2)
%interptyp=1;
%interptyp=input('Interpolation Type [1 = MR Voxels, 2 = Mesh Locations] (default = 1): ');
%if isempty(interptyp)
    interptyp=1;
%end

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

d1=dir('../*.nod');
% remove '.hom.nod', regSP.nod and reg.nod entries (just take shortest)
%remove '.hom.nod', regSP.nod and reg.nod entries (just take shortest)
ifle=0;
for ii=1:length(d1)
    if (length(d1(ii).name)>7)&&(strcmp(d1(ii).name(end-7:end),'.hom.nod'))
        %ignore this file
    elseif(length(d1(ii).name)>7)&&(strcmp(d1(ii).name(end-7:end),'.reg.nod'))
        % ignore this file
    elseif((length(d1(ii).name)>9)&&(strcmp(d1(ii).name(end-9:end),'.regSP.nod')))
        %ignore this file
    else
        ifle=ifle+1;
        dnod(ifle).name=d1(ii).name;          
    end
end


% ifle=0;
% for ii=1:length(d1)
%     if(length(d1(ii).name)>9)
%         if( (~strcmp(d1(ii).name(end-7:end),'.hom.nod')) &&(~strcmp(d1(ii).name(end-7:end),'.reg.nod'))&&(~strcmp(d1(ii).name(end-9:end),'.regSP.nod')))
%             ifle=ifle+1;
%             dnod(ifle).name=d1(ii).name;          
%         end
%     end
% end

%displacement mesh
for ii=1:length(dnod)
    disp(['File ' int2str(ii) ' :: ' dnod(ii).name])
end
clear n
if(length(dnod)==1)
    n=1;
else
    n=input('Number of file to use (default = last)  >> ');
end

if(isempty(n))
    n=length(dnod);
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
    
    if(size(junk,2)>4)
        mesh(ii).nodsense=junk(:,5);
    else
        mesh(ii).nodsense=ones(size(junk,1),1);
    end
    
    if(size(junk,2)>5)
        mesh(ii).nodreg=junk(:,6);
    else
        mesh(ii).nodreg=zeros(size(junk,1),1);
    end
    
    
    if(exist([d(ii).name(1:end-3) 'elm'],'file'))
        junk=load([d(ii).name(1:end-3) 'elm']);
        mesh(ii).in=junk(:,2:9);
        mesh(ii).ne=length(junk);
        
        if(size(junk,2)>28)
            mesh(ii).elmsense=junk(:,29);
        else
            mesh(ii).elmsense=ones(size(junk,1),1);
        end
        
        if(size(junk,2)>29)
            mesh(ii).elmreg=junk(:,30);
        else
            mesh(ii).elmreg=zeros(size(junk,1),1);
        end
        
        
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
    
    if~(strcmp(msk(end-3:end),'.mat'))
        msk=[msk '.mat'];
    end
    disp(['Loading Mask and MRE3DMotionData.mat File:  ' '../' msk])
    mskf=['../../../' msk];
    if exist(mskf,'file')
        load(mskf);
    else
        mskf=input('Default Mask File Not Found, Please Enter Path and Name:  >>','s');
        load(mskf);
    end
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

%% Plot properties

displayinterp=1; % Factor to interpolate images for display (set to 2 to double the size, sometimes the images are too small on screen to see.)

if ~meshmag
clear intstack
for jj=1:size(MagIm,3)
    if(displayinterp==1)
        intstack(:,:,jj)=MagIm(:,:,jj);
    else
        intstack(:,:,jj)=interp2(MagIm(:,:,jj),displayinterp);
    end
end
%save('intstack.mat','intstack');
if(showimg)
    montagestack(intstack);
    title('MR Magnitude image')
end
else
    MagImStack=interp3(nodY27,nodX27',nodZ27,MagIm,Ymr,Xmr',Zmr,'linear',0).*mask;
    clear intstack
    for jj=1:size(MagImStack,3)
        if(displayinterp==1)
            intstack(:,:,jj)=MagImStack(:,:,jj);
        else
            intstack(:,:,jj)=interp2(MagImStack(:,:,jj),displayinterp);
        end
    end
    %save('intstack.mat','intstack');
    if(showimg)
        montagestack(intstack);
        title('MR Magnitude image')
    end
end

d=dir([mtrstm 'RE..' itrstr '.prop.*.mtr']);
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
        prop(ii).revals=load([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
        prop(ii).imvals=load([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
    else
        prop(ii).revals=load([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
        prop(ii).imvals=load([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
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
    
    for jj=1:prop(ii).npr
        if(basis(ii)==1)
            i1=mesh(rmsh).ind(jj,1);
            i2=mesh(rmsh).ind(jj,2);
            i3=mesh(rmsh).ind(jj,3);
            sns=mesh(rmsh).nodsense(jj);
        elseif(basis(ii)==2)
            i1=mesh(rmsh).indelm(jj,1);
            i2=mesh(rmsh).indelm(jj,2);
            i3=mesh(rmsh).indelm(jj,3);
            sns=mesh(rmsh).elmsense(jj);
        end
        if(sns==0) % Materialprop point is not updated
            prop(ii).restack(i1,i2,i3)=NaN;
            if mfpropind(ii)
                prop(ii).restack(i1,i2,i3,2)=NaN;
            end
        else
            prop(ii).restack(i1,i2,i3)=prop(ii).revals(jj,2);
            if mfpropind(ii)
                prop(ii).restack(i1,i2,i3,2)=prop(ii).revals(jj,3);
            end
        end
        
    end
    for jj=1:prop(ii).npi
        if(basis(ii)==1)
            i1=mesh(imsh).ind(jj,1);
            i2=mesh(imsh).ind(jj,2);
            i3=mesh(imsh).ind(jj,3);
            sns=mesh(imsh).nodsense(jj);
        elseif(basis(ii)==2)
            i1=mesh(imsh).indelm(jj,1);
            i2=mesh(imsh).indelm(jj,2);
            i3=mesh(imsh).indelm(jj,3);
            sns=mesh(imsh).elmsense(jj);
        end
        if(sns==0) % Materialprop point is not updated
            prop(ii).imstack(i1,i2,i3)=NaN;
            if mfpropind(ii)
                prop(ii).imstack(i1,i2,i3,2)=NaN;
            end
        else
            prop(ii).imstack(i1,i2,i3)=prop(ii).imvals(jj,2);
            if mfpropind(ii)
                prop(ii).imstack(i1,i2,i3,2)=prop(ii).imvals(jj,3);
            end
        end
        
        
    end
    
    %SplineInterp3D_withnans(xin,yin,zin,vin,xout,yout,zout,btyp,contspline)
    % Interpolate to Chosen Resolution (MR resolution, interptyp=1, or mesh resolution, interptyp=2)
    if(basis(ii)==1)
        %prop(ii).reMRstack=interp3(mesh(rmsh).nodY,mesh(rmsh).nodX',mesh(rmsh).nodZ,prop(ii).restack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        %prop(ii).imMRstack=interp3(mesh(imsh).nodY,mesh(imsh).nodX',mesh(imsh).nodZ,prop(ii).imstack(:,:,:,1),Ymr,Xmr',Zmr,'linear',0).*mask;
        prop(ii).reMRstack=LinearInterp3D_withnans(mesh(rmsh).nodX,mesh(rmsh).nodY,mesh(rmsh).nodZ,prop(ii).restack(:,:,:,1),Xmr,Ymr,Zmr);
        prop(ii).imMRstack=LinearInterp3D_withnans(mesh(imsh).nodX,mesh(imsh).nodY,mesh(imsh).nodZ,prop(ii).imstack(:,:,:,1),Xmr,Ymr,Zmr); 
        
        if mfpropind(ii)
            prop(ii).reMRstack(:,:,:,2)=LinearInterp3D_withnans(mesh(rmsh).nodX,mesh(rmsh).nodY,mesh(rmsh).nodZ,prop(ii).restack(:,:,:,2),Xmr,Ymr,Zmr);
            prop(ii).imMRstack(:,:,:,2)=LinearInterp3D_withnans(mesh(imsh).nodX,mesh(imsh).nodY,mesh(imsh).nodZ,prop(ii).imstack(:,:,:,2),Xmr,Ymr,Zmr);
        end
    elseif(basis(ii)==2)
        prop(ii).reMRstack=LinearInterp3D_withnans(mesh(rmsh).elmX,mesh(rmsh).elmY,mesh(rmsh).elmZ,prop(ii).restack(:,:,:,1),Xmr,Ymr,Zmr,2);
        prop(ii).imMRstack=LinearInterp3D_withnans(mesh(imsh).elmX,mesh(imsh).elmY,mesh(imsh).elmZ,prop(ii).imstack(:,:,:,1),Xmr,Ymr,Zmr,2);
        if mfpropind(ii)
            prop(ii).reMRstack(:,:,:,2)=LinearInterp3D_withnans(mesh(rmsh).elmX,mesh(rmsh).elmY,mesh(rmsh).elmZ,prop(ii).restack(:,:,:,2),Xmr,Ymr,Zmr);
            prop(ii).imMRstack(:,:,:,2)=LinearInterp3D_withnans(mesh(imsh).elmX,mesh(imsh).elmY,mesh(imsh).elmZ,prop(ii).imstack(:,:,:,2),Xmr,Ymr,Zmr);
        end
    end
    
    prop(ii).reMRstack(mask==0)=NaN;
    prop(ii).imMRstack(mask==0)=NaN;
    if mfpropind(ii)
        junkstack=prop(ii).reMRstack(:,:,:,2);
        junkstack(mask==0)=NaN;
        prop(ii).reMRstack(:,:,:,2)=junkstack;
        junkstack=prop(ii).imMRstack(:,:,:,2);
        junkstack(mask==0)=NaN;
        prop(ii).imMRstack(:,:,:,2)=junkstack;
    end
    
    
    if(rcnind(2*(ii-1)+1))
        clear intstack
        for jj=1:size(prop(ii).reMRstack(:,:,:,1),3)
           if(displayinterp==1)
               intstack(:,:,jj)=prop(ii).reMRstack(:,:,jj,1);
           else
               intstack(:,:,jj)=interp2(prop(ii).reMRstack(:,:,jj,1),displayinterp);
           end
        end
        if(showimg)
            montagestack(intstack);
            colorbar
            title(['Real Prop ' int2str(ii) '; res = ' num2str(mesh(rmsh).res)])
        end
        if mfpropind(ii)
            clear intstack
            for jj=1:size(prop(ii).reMRstack(:,:,:,2),3)
                if(displayinterp==1)
                    intstack(:,:,jj)=prop(ii).reMRstack(:,:,jj,2);
                else
                    intstack(:,:,jj)=interp2(prop(ii).reMRstack(:,:,jj,2),displayinterp);
                end
            end
            if(showimg)
                montagestack(intstack); 
                colorbar
                title(['Real Prop ' int2str(ii) ' \alpha component; res = ' num2str(mesh(rmsh).res)])
            end
        end
    end
    
    if(rcnind(2*(ii-1)+2))
        clear intstack
        for jj=1:size(prop(ii).imMRstack(:,:,:,1),3)
            if(displayinterp==1)
                intstack(:,:,jj)=prop(ii).imMRstack(:,:,jj,1);
            else
                intstack(:,:,jj)=interp2(prop(ii).imMRstack(:,:,jj,1),displayinterp);
            end
        end
        if(showimg);
            montagestack(intstack); 
            colorbar
            title(['Imag Prop ' int2str(ii) '; res = ' num2str(mesh(imsh).res)])
        end
        if mfpropind(ii)
            clear intstack
            for jj=1:size(prop(ii).imMRstack(:,:,:,2),3)
                if(displayinterp==1)
                    intstack(:,:,jj)=prop(ii).imMRstack(:,:,jj,2);
                else
                    intstack(:,:,jj)=interp2(prop(ii).imMRstack(:,:,jj,2),displayinterp);
                end
            end
            if(showimg);
                montagestack(intstack);
                colorbar
                title(['Imag Prop ' int2str(ii) ' \alpha component; res = ' num2str(mesh(imsh).res)])
            end
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

if(mfind) 
    dfreq=input('Frequency for Damping Property Display: ');
end

if ~mfpropind(1)
    murat=prop(1).imMRstack(:,:,:,1)./prop(1).reMRstack(:,:,:,1);
else
    murat=(prop(1).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).imMRstack(:,:,:,2)))./(prop(1).reMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).reMRstack(:,:,:,2)));
    mure=(prop(1).reMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).reMRstack(:,:,:,2)));
    muim=(prop(1).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(1).imMRstack(:,:,:,2)));
end
if ~mfpropind(2)
    if rcnind(4)
        rhorat=prop(2).imMRstack(:,:,:,1)./prop(2).reMRstack(:,:,:,1);
    else
        rhorat=0;
    end
else
    rhorat=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)))./prop(2).reMRstack(:,:,:,1);
    rhoim=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)));
end
DR=0.5.*(murat-rhorat);
RC=0.5*murat./DR;

if mfpropind(1)
    clear intstack
    for jj=1:size(mure,3)
        if(displayinterp==1)
            intstack(:,:,jj)=mure(:,:,jj);
        else
            intstack(:,:,jj)=interp2(mure(:,:,jj),displayinterp);
        end
    end
    if(showimg);
        montagestack(intstack);
        title(['\mu_R at ',num2str(dfreq),' Hz']);
        colorbar
    end
    
    clear intstack
    for jj=1:size(muim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=muim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(muim(:,:,jj),displayinterp);
        end
    end
    if(showimg); 
        montagestack(intstack)
        title(['\mu_I at ',num2str(dfreq),' Hz']);
        colorbar
    end
end

if mfpropind(2)
    clear intstack
    for jj=1:size(rhoim,3)
        if(displayinterp==1)
            intstack(:,:,jj)=rhoim(:,:,jj);
        else
            intstack(:,:,jj)=interp2(rhoim(:,:,jj),displayinterp);
        end
    end
    if(showimg); 
        montagestack(intstack)
        title(['\rho_I at ',num2str(dfreq),' Hz']);
        colorbar
    end
end

clear intstack
for jj=1:size(DR,3)
    if(displayinterp==1)
        intstack(:,:,jj)=DR(:,:,jj);
    else
        intstack(:,:,jj)=interp2(DR(:,:,jj),displayinterp);
    end
end
if(showimg); 
    montagestack(intstack);
    if ~mfpropind(1)
        title('Damping Ratio')
    else
        title(['Damping Ratio at ',num2str(dfreq),' Hz']);
    end
    colorbar
end

clear intstack
for jj=1:size(RC,3)
    if(displayinterp==1)
        intstack(:,:,jj)=RC(:,:,jj);
    else
        intstack(:,:,jj)=interp2(RC(:,:,jj),displayinterp);
    end
end
if(showimg)
    montagestack(intstack)
    if ~mfpropind(1)
        title('Rayleigh Composition')
    else
        title(['Rayleigh Composition at ',num2str(dfreq),' Hz']);
    end
    colorbar
end

%% Read and process .invparam file
invparamf=[nodstm 'invparam'];
if(exist(invparamf,'file'))
    disp('Inversion Parameters stored in invparam structure')
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
          meshf=input('Default Mesh Input File Not Found, Please Enter Path and Name:  >>','s');
          load(meshf);
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

% Save figures
if(showimg)
    figsave=input('Format to save figures (blank=no save) : >> ','s');


    if(isempty(figsave))
        disp('- No figures saved')
    else
        figHandles = findobj('Type','figure');
        disp('Saving Figures')
        for ii=1:length(figHandles)
            saveas(figHandles(ii),[mtrstm '.' itrstr 'reconimg' int2str(figHandles(ii)) '.' figsave],figsave)
        end
    end
end
    
    



    
    








