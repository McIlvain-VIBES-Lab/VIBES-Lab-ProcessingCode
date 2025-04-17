% MRE_plotv9_mask_Linear_Function: Converts MREv8 property files into a 
%    3D matlab stack on the same points as MRE data.
%    Works independently of the type of displacement mesh
%    Should be run from the 'inv' folder. Assumes a specific folder
%    structure:
%    Basename/meshtyp/meshname/inv
%    Directory Basename contains MRE_3DMotionData.mat and HeaderData.mat file
%    Directory  meshtyp is hex or tet
%    Directory meshname is the name of the mesh
%    Directory inv contains the appropriate meshind, mtr and property mesh
%    flies (.mtrmesh.##.nod and mtrmesh.##.elm)


% Outputs results as a mat file
% Inputs: 
%   Results are interpolated back to the original MRE data resolution.
%       - Records reconstuction and mesh data in the same structure as the 
%         results.
%       - Uses data saved from the meshing step to interpolate material
%         properties back onto the original MR voxels. 
%close all;clear all;clc

function MRE_plotv9_mask_Linear_Function(nodf,mtrf)


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

basisf=[nodstm 'basis'];
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
    basis=dlmread(basisf);
else
    basis=[1 1 1];
end
if (sum(basis==2)>0)
    elmbasis=true;
else
    elmbasis=false;
end

d=dir([nodstm 'mtrmesh*nod']);
nmesh=length(d);
for ii=1:nmesh
    junk=dlmread(d(ii).name);
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
        junk=dlmread([d(ii).name(1:end-3) 'elm']);
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
      %disp('XXX NEED TO CONSIDER BASIS HERE')
      if(jj==1) 
          mesh(ii).nodX=srt;          
      elseif(jj==2)
          mesh(ii).nodY=srt;          
      elseif(jj==3)          
          mesh(ii).nodZ=srt;          
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

% These are the mesh resolution, which may be interpolated from the
% original data

% Load imginfo file from one directory back
imgf=dir('../*imginfo');
imgf=imgf(1); % In case there are 2 meshes in the folder
if(isempty(imgf))
    [imgf(1).name, pathname] = uigetfile('../*imginfo', 'Select appropriate imginfo or MRE_3DMotionData.mat file for interpolation');
    if(strcmp(imgf(1).name,'MRE_3DMotionData.mat'))
        load([pathname 'MRE_3DMotionData.mat']);
        if(exist([pathname 'HeaderData.mat'],'file'))
            load([pathname 'HeaderData.mat']);
            mridim=size(MagIm);
            voxsize_mm=DirIndex(4,1:3);
            RHcoord=DirIndex(3,3)~=-1;
        elseif(exist([pathname 'MREMotionInfo.mat'],'file'))
            load([pathname 'MREMotionInfo.mat']);
            mridim=[MREMotionInfo.nX MREMotionInfo.nY MREMotionInfo.nS];
            voxsize_mm=[MREMotionInfo.pixelspacing_x MREMotionInfo.pixelspacing_y MREMotionInfo.slicethickness+MREMotionInfo.slicegap];
            RHcoord=true;
        end
        
        % Write imginfo file
        
        junk=dir('../*bnod');
        fid=fopen(['../' junk(1).name(1:end-4) '.imginfo'],'w');
          fprintf(fid,'Size of MRI stack \n');
          fprintf(fid,'%7i %7i %7i \n',mridim);
          fprintf(fid,'Voxel dimensions (mm) \n');
          fprintf(fid,'%15.6e %15.6e %15.6e \n',voxsize_mm);
          fprintf(fid,'Right hand coordinate system indicator \n');
          fprintf(fid,'%i',RHcoord);
        fclose(fid);
        
        imgf.name=[junk(1).name(1:end-5) '.imginfo'];
    end
end
fid=fopen(['../' imgf.name]);
tline=fgetl(fid);
tline=fgetl(fid);
mridim=str2num(tline);
tline=fgetl(fid);
tline=fgetl(fid);
voxsize_mm=str2num(tline);
tline=fgetl(fid);
tline=fgetl(fid);
RHcoord=str2num(tline);
fclose(fid)


disp(['Loading HeaderData.mat and MRE3DMotionData.mat File:'])
% Get MagIm Values from MRE Data
if exist(['../../../MRE_3DMotionData.mat'],'file')
    load(['../../../MRE_3DMotionData.mat'])    
elseif exist(['../../MRE_3DMotionData.mat'],'file')
    load(['../../MRE_3DMotionData.mat'])    
elseif exist(['../../../NLI_MRE_Input.mat'],'file')
    load(['../../../NLI_MRE_Input.mat']);
    MagIm=AnatomicalMRI;
elseif exist(['../../NLI_MRE_Input.mat'],'file')
    load(['../../NLI_MRE_Input.mat']) ;   
    MagIm=AnatomicalMRI;
else
    warning('cant find MRE_3DMotionData.mat or NLI_MRE_Input.mat - no MagIm info available')
    MagIm=zeros(mridim);
end
sz=size(MagIm);
Xmr=(0:sz(1)-1).*voxsize_mm(1)/1000;  % May be wrong for anisotropic voxels - DirIndex transformation not considered
Ymr=(0:sz(2)-1).*voxsize_mm(2)/1000;
Zmr=(0:sz(3)-1).*voxsize_mm(3)/1000;

if(~RHcoord) % Undo the reversal of the z direction for non RH coord systems
    disp('Undoing Z direction revesal for non RH coord system');
    for ii=1:length(mesh)
        mesh(ii).nodZ=-mesh(ii).nodZ;
    end
end

%% Plot properties
d=dir([mtrstm 'RE..' itrstr '.prop.*.mtr']);
numprop=length(d);
mfpropind=false(numprop,1);
meshind=zeros(2,numprop);
for ii=1:numprop
    if mfind
        if d(ii).name(end-5:end-4)=='mf'
            mfpropind(ii)=true;
        end
    end
    rmsh=dlmread([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.meshind']);
    imsh=dlmread([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.meshind']);
    meshind(1,ii)=rmsh;
    meshind(2,ii)=imsh;
    if ~mfpropind(ii)
        prop(ii).revals=dlmread([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
        prop(ii).imvals=dlmread([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mtr']);
    else
        prop(ii).revals=dlmread([mtrstm 'RE..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
        prop(ii).imvals=dlmread([mtrstm 'IM..' itrstr '.prop.' sprintf('%2.2i',ii) '.mf.mtr']);
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
    
    
end
%%


%%



% Poisson Ratio


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
                  disp(['invparam.' strtrim(vname) ' = [' values '];'])
                  eval(['invparam.' strtrim(vname) ' = [' values '];'])   
              end
          else
              endoffile=true;
          end
          
      end
      if(exist(['../' invparam.nodefile(1:end-7) 'meshinput.mat'],'file'))
          meshinputs=load(['../' invparam.nodefile(1:end-7) 'meshinput.mat']);
      else
          meshinputs='notrecorded';
      end
    fclose(fid); 
else
    invparam='Data Not Present - MREv7.04 and above required for recording of inversion parameters';
end

if(exist(invparamf,'file')) % Info about the material model is present
    modeltyp=invparam.modeltype;
    if(isfield(invparam,'nod_per_disp_elm'))
        npe=invparam.nod_per_disp_elm;
    else
        npe=27;
    end
else % Assume it is model = 1 (iso incompressible, all versions where invparamf is not present can only handle isocompressible
    modeltyp=1;
    npe=27;
end

% Build damping ratio images. 
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
    if ((modeltyp==1)||(modeltyp==3))
        rhorat=prop(2).imMRstack(:,:,:,1)./prop(2).reMRstack(:,:,:,1);
    elseif((modeltyp==5)||(modeltyp==6))
        rhorat=prop(4).imMRstack(:,:,:,1)./prop(4).reMRstack(:,:,:,1);
    else
        rhorat=0;
    end
else
    if ((modeltyp==1)||(modeltyp==3))
      rhorat=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)))./prop(4).reMRstack(:,:,:,1);
      rhoim=(prop(2).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(2).imMRstack(:,:,:,2)));
    elseif((modeltyp==5)||(modeltyp==6))
      rhorat=(prop(4).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(4).imMRstack(:,:,:,2)))./prop(4).reMRstack(:,:,:,1);
      rhoim=(prop(4).imMRstack(:,:,:,1).*((2*pi*dfreq).^prop(4).imMRstack(:,:,:,2)));
    else
        rhorat=0;
    end
end
DR=0.5.*(murat-rhorat);
RC=0.5*murat./DR;

%% Save properties

if(modeltyp==1) % iso incompressible 
    RealShear=prop(1).reMRstack;
    ImagShear=prop(1).imMRstack;
    RealDens=prop(2).reMRstack;
    ImagDens=prop(2).imMRstack;
    RealBulk=prop(3).reMRstack;
    ImagBulk=prop(3).imMRstack;    
    Poisson=(3*(prop(3).reMRstack+1i*prop(3).imMRstack)-2*(prop(1).reMRstack+1i*prop(1).imMRstack)) ./ (2*(3*(prop(3).reMRstack+1i*prop(3).imMRstack)+(prop(1).reMRstack+1i*prop(1).imMRstack)));  
    RealPoisson=real(Poisson).*(RealShear>0);
    ImagPoisson=imag(Poisson).*(RealShear>0);
    ShearModMag=abs(prop(1).reMRstack + 1i*prop(1).imMRstack);
    ShearModPhaseAngle=angle(prop(1).reMRstack + 1i*prop(1).imMRstack);
    ShearStiffness= 2 * ShearModMag.^2 ./ (prop(1).reMRstack + ShearModMag);
elseif(modeltyp==3) % iso incompressible 
    RealShear=prop(1).reMRstack;
    ImagShear=prop(1).imMRstack;
    RealDens=prop(2).reMRstack;
    ImagDens=prop(2).imMRstack;
    RealLambda=prop(3).reMRstack;
    ImagLambda=prop(3).imMRstack;    
    Poisson=(prop(3).reMRstack+1i*prop(3).imMRstack) ./ (2*((prop(3).reMRstack+1i*prop(3).imMRstack) + (prop(1).reMRstack+1i*prop(1).imMRstack))); 
    RealPoisson=real(Poisson).*(RealShear>0);
    ImagPoisson=imag(Poisson).*(RealShear>0);
    ShearModMag=abs(prop(1).reMRstack + 1i*prop(1).imMRstack);
    ShearModPhaseAngle=angle(prop(1).reMRstack + 1i*prop(1).imMRstack);
    ShearStiffness= 2 * ShearModMag.^2 ./ (prop(1).reMRstack + ShearModMag);
elseif(modeltyp==4) % Poroelastic 
    RealShear=prop(1).reMRstack;
    ImagShear=prop(1).imMRstack;
    RealLambda=prop(2).reMRstack;
    ImagLambda=prop(2).imMRstack;
    RealHC=prop(3).reMRstack;
    ImagHC=prop(3).imMRstack;
    RealPorosity=prop(4).reMRstack;
    ImagPorosity=prop(4).imMRstack;    
    Poisson=(prop(2).reMRstack+1i*prop(2).imMRstack) ./ (2*((prop(2).reMRstack+1i*prop(2).imMRstack) + (prop(1).reMRstack+1i*prop(1).imMRstack))); 
    RealPoisson=real(Poisson).*(RealShear>0);
    ImagPoisson=imag(Poisson).*(RealShear>0);
elseif(modeltyp==5) % NITI 3 moduli
    Realmu=prop(1).reMRstack;
    Imagmu=prop(1).imMRstack;
    Realmus=prop(2).reMRstack;
    Imagmus=prop(2).imMRstack;
    Realmut=prop(3).reMRstack;
    Imagmut=prop(3).imMRstack;
    RealDens=prop(4).reMRstack;
    ImagDens=prop(4).imMRstack;    
    RealBulk=prop(5).reMRstack;
    ImagBulk=prop(5).imMRstack;    
elseif(modeltyp==6) % NITI mu phi zeta 
    Realmu=prop(1).reMRstack;
    Imagmu=prop(1).imMRstack;
    Realphi=prop(2).reMRstack;
    Imagphi=prop(2).imMRstack;
    Realzeta=prop(3).reMRstack;
    Imagzeta=prop(3).imMRstack;
    RealDens=prop(4).reMRstack;
    ImagDens=prop(4).imMRstack;    
    RealBulk=prop(5).reMRstack;
    ImagBulk=prop(5).imMRstack;  
    
end






for ii=1:nmesh
  meshres(ii,:)=mesh(ii).res;
end

outmatf=[mtrstm '.' itrstr '.ReconProps.mat'];
recondir=pwd;
save(outmatf,'MagIm','meshind','meshres','recondir');
if(exist(invparamf,'file')) % Info about the inversion params is present
    save(outmatf,'meshinputs','invparam','-append')
end
if(modeltyp==1) % iso incompressible 
    save(outmatf,'RealShear','ImagShear','RealDens','ImagDens','RealBulk','ImagBulk','DR','RC','ShearStiffness','ShearModMag','ShearModPhaseAngle','RealPoisson','ImagPoisson','-append')
elseif(modeltyp==3) % iso compressible 
    save(outmatf,'RealShear','ImagShear','RealDens','ImagDens','RealLambda','ImagLambda','DR','RC','ShearStiffness','ShearModMag','ShearModPhaseAngle','RealPoisson','ImagPoisson','-append')
elseif(modeltyp==4) % Poroelastic
    save(outmatf,'RealShear','ImagShear','RealLambda','ImagLambda','RealHC','ImagHC','RealPorosity','ImagPorosity','RealPoisson','ImagPoisson','-append')
elseif(modeltyp==5) % NITI 3 moduli
    save(outmatf,'Realmu','Imagmu','Realmus','Imagmus','Realmut','Imagmut','RealDens','ImagDens','RealBulk','ImagBulk','-append')
elseif(modeltyp==6) % NITI standard
    save(outmatf,'Realmu','Imagmu','Realphi','Imagphi','Realzeta','Imagzeta','RealDens','ImagDens','RealBulk','ImagBulk','DR','RC','-append')

end


end
    
    



    
    








