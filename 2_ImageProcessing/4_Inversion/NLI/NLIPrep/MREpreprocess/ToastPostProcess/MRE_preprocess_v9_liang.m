function MRE_preprocess_v9_liang(default,SubjectName)
% This is a modified version. 
% The following lines were modified, starting from the line 1234.
% if(strcmpi(meshtype,'hex')) % Output Hexahedral runfiles
%     output_runfiles_hex_v9_liang(outstm,vox,muest,rhoest,DR,inpath,nodhmgf,elmhmgf,bcoutf,regoutf,freqHz,dspoutf,outpath,znedgelength,znovlp,noreg,noDTI,SubjectName);
% elseif(strcmpi(meshtype,'tet'))
%     output_runfiles_tet_v9_liang(outstm,vox,muest,rhoest,DR,inpath,nodhmgf,elmhmgf,bcoutf,regoutf,pbcf2,freqHz,dspoutf,outpath,znedgelength,znovlp,noreg,noDTI,SubjectName);
% end


% % MRE_preprocess Generates required files for MREv8 from MRE data. 
%   Interpolates data to an approximate number of nodes per wavlength,
%   based on prior estiamtes of the stiffness. Builds Hex27 mesh manually,
%   i.e. does not require a 'template mesh'.
%   If called using MRIhexmesh('default'), uses the default values without
%   prompting for inputs. 
%   INPUTS:
%    default(optional): Set to 'default' to disable propmts for inputs. Any
%                       other value, or if not supplied will propmt for
%                       inputs. Using default values allows automated
%                       meshing.
%
%   IMPORTANT FEATURES
%
%   - DISPLACEMENT DATA IS SCALED!! I chose to scale everything so that the
%     average displacement amplitude is the same for all datasets
%     interpolated with this code. 
%   - Interpolation is an option. Cubic spline interpoaltion is used, a
%     custom function is used which does not use values outside of the mask
%     for the interpolated values.
%   - Different meshing strategies are possible: (MRE-Zonev7.05 now takes
%         zone dimensions rather than dividing the extents into an integer
%         number of zones.)
%      meshstrat=1 uses one FE node for each MR voxel, as has been done in
%                  the past.
%      meshstrat=2 allows specification of the number of FE nodes per
%                  wavelength, based on an estiamte of the shear modulus. A
%                  warning will be produced if the data is being
%                  interpolated to a resolutionlower than the data.
%      meshstrat=3 allows specification of a target resolution.
%      Both meshstrat=2 and meshstrat=3 make slight modifications to the
%      resolution to fit an integer number of 27 node elements across the
%      domain defined by the extents of the mask, to avoid wasting planes
%      of data near the edges.
%   - Different zone sizing strategies are possible: 
%      zonestrat=1 Matches a number of FE nodes per zone, as has been done 
%                  in the past.
%      zonestrat=2 matches a number of wavelengths per zone, based on an 
%                  estiamte of the shear modulus 
%      zonestrat=3 define a value for the zone size
%      Both of these strategies are only approximate, the way that the
%      zoning proscess works (dividing the mesh extents into an integer
%      number of zones) makes it difficult to exactly match specified
%      values.
%   - Input data which may be useful when analyzing results or looking at
%     a reconstruction later on is saved as outstm.InterpLocations.mat,
%     and outstm.InterpData.mat outstm.meshinput.mat. The idea of saving
%     this data is to attach it to reconstruction results so it is clear
%     what parameters were used to produce it. 
%   
disp('MRE-Zone v9 Data Converter')
par=false; % If you have the matlab parallel processing toolbox set this to true, otherwise false.
nlabs=7; % Number of labs for matlab to use (max 7)
if(par)
    disp(['Using Parallelized version with ' int2str(nlabs) ' labs'])
else
    disp('Using non-Parallelized version - modify value of ''par'' to change')
end

if(nargin==0) % Default value is to prompt for inputs.
    default='no';
end
usedef=strcmp(default,'default');
if(usedef) 
    disp('Using Default values, no input prompts')
end

%% Get Inputs

% MRE data files:
mridef='NLI_MRE_Input.mat';
if(~usedef) 
    MRIfile=input(['NLI-MRE input file name? (default is ' mridef '):  >>'],'s');
end
if ~exist('MRIfile','var')||isempty(MRIfile)
    MRIfile=mridef;
end

mskdef='Mask.mat';
if(~usedef) 
    msk=input(['Name of mask to use for mesh (Default ' mskdef '):  >>'],'s');
end 
if ~exist('msk','var')||isempty(msk)
    msk=mskdef;
end

if(exist('SPRregs.mat','file'))
    regdef='SPRregs.mat';
else
    regdef='none';
end
if(~usedef) 
    regf=input(['Name of Segmented Stack File (Default ' regdef '):  >>'],'s');
end 
if ~exist('regf','var')||isempty(regf)
    regf=regdef;
end

if(exist('DTI.mat','file'))
    DTIdef='DTI.mat';
else
    DTIdef='none';
end
if(~usedef) 
    DTIf=input(['Name of DTI Fiber Directions File (Default ' DTIdef '):  >>'],'s');
end 
if ~exist('DTIf','var')||isempty(DTIf)
    DTIf=DTIdef;
end



% Mesh type (27 node hexahedral or 4 node tetrahedral)
meshtypedef='hex';
if(~usedef) 
    disp('Mesh Element Type')
    meshtype=input(['Element type <hex,tet>, (Default ' meshtypedef ') >>'],'s');
end
if ~exist('meshtype','var')||isempty(meshtype)
    meshtype=meshtypedef;close all
    
end

%% Load and check the input data

load(MRIfile);

% Check data format:
if(~exist('Motion2Image','var')) % Backward compatibility:
    if(exist('A','var'))
        warning('NLI MRE data in old format, running backward compatibility script')
        [Ur,Ui,Motion2Image,RHcoord,voxsize_mm,freqHz,AnatomicalMRI]=MRE_MotionData_convert(MRIfile);
        clear A P MagIm
    else
        error('NLI MRE input in unrecognized format')
    end
end


if(strcmpi(meshtype,'hex'))
    %Meshing Approach: one node per voxel, or interpolate to give nodes per wavelength
    meshstratdef=1;
    if(~usedef) 
        disp(['Meshing Strategies: 1=one node per MR voxel, 2=Target #nodes per wavelength, 3=Target Resolution'])
        meshstrat=input(['Meshing Strategy, (Default ' int2str(meshstratdef) ') >>']);
    end
    if ~exist('meshstrat','var')||isempty(meshstrat)
        meshstrat=meshstratdef;
    end
elseif(strcmpi(meshtype,'tet'))
    %Meshing Approach: one node per voxel, or interpolate to give nodes per wavelength
    tetmeshresdef=mean(voxsize_mm);
    if(~usedef)         
        tetmeshres=input(['Tetrahedral Mesh Resolution (mm), (Default ' num2str(tetmeshresdef) ') >>']);
    end
    if ~exist('tetmeshres','var')||isempty(tetmeshres)
        tetmeshres=tetmeshresdef;
    end
    
    tetcustomdef=false;
    if(~usedef)
        tetcustom=input('Do you want to customize your tetrahedral meshing criteria? 1=yes, (Default: no)>>:');
    end
    if ~exist('tetcustom','var')||isempty(tetcustom)
        tetcustom=tetcustomdef;
    end
    
else
    error(['Meshtype :: ' meshtype ' not recognized'])
end


%Zone Sizing Approach: nodes per zone or wavelngths per zone
zonestratdef=3;
if(~usedef) 
    disp(['Zone sizing Strategies: 1=fit nodes per zone, 2=fit wavelengths per zone, 3 = fixed size'])
    zonestrat=input(['Zone sizing Strategy, (Default ' int2str(zonestratdef) ') >>']);
end
if ~exist('zonestrat','var')||isempty(zonestrat)
    zonestrat=zonestratdef;
end

mudef=3300;
%mudef=20000;
% Estimated properties.    
if(~usedef) 
    muest=input(['Estimate of Shear Modulus (Default ' num2str(mudef) ') >>']);
end
if ~exist('muest','var')||(isempty(muest))
    muest=mudef;
end    

DRdef=0.18;
if(~usedef) 
    DR=input(['Estimate of Damping ratio (Default ' num2str(DRdef) ') >>']);
end
if ~exist('DR','var')||(isempty(DR))
    DR=DRdef;
end    




load(msk);

erodedef=0;
if(~usedef) 
    erodemaskn=input(['Number of times to erode mask ?? default: ' int2str(erodedef) ' ) >> ']);
end
if ~exist('erodemaskn','var')||(isempty(erodemaskn))
    erodemaskn=erodedef;
end   

cleandef='n';
if(~usedef) 
    cleanmaskflg=input(['Clean hanging mask voxels ?? <y n> default: ' cleandef ' ) >> '],'s');
end
if ~exist('cleanmaskflag','var')||(isempty(cleanmaskflg))
    cleanmaskflg=cleandef;
end    

if(strcmpi(meshtype,'tet'))
    smoothdef='n';
    if(~usedef) 
        smoothmesh=input(['Smooth tet mesh surface ?? <y n> default: ' smoothdef ' ) >> '],'s');
    end
    if ~exist('smoothmesh','var')||(isempty(smoothmesh))
        smoothmesh=smoothdef;
    end    
    
    if(strcmp(smoothmesh,'y'))
       smoothfacdef=3;
        if(~usedef) 
            smoothfac=input(['smoothing paramter ?? <y n> default: ' int2str(smoothfacdef) ' ) >> ']);
        end
        if ~exist('smoothfac','var')||(isempty(smoothfac))
            smoothfac=smoothfacdef;
        end 
    end
end

if(erodemaskn>0)
    disp(['Eroding mask in-plane ' int2str(erodemaskn) ' times'])
    for ii=1:erodemaskn
        for jj=1:size(mask,3)
            mask(:,:,jj)=imerode(mask(:,:,jj),ones(3,3));
        end
    end
end

if(strcmpi(cleanmaskflg,'y'))
    mask=cleanmask(mask);
end
    
    





if(ndims(AnatomicalMRI)==4);
    AnatomicalMRI=mean(AnatomicalMRI,4);
end

% Note that for non-axisymmetric data vox should be the sizes in the  [ (:,1,1) (1,:,1) (1,1,:) ] directions
vox=voxsize_mm.*1e-3;

mridim=size(mask);


rhoest=1000;

wl=1/freqHz*sqrt(muest/rhoest); % Undamped elastic wavelength, could include full VE shear wavelength if neccesary 

if(strcmpi(meshtype,'hex'))
    if(meshstrat==2)
        % Desired nodes per wavelength.
        npwdef=12;
        if(~usedef) 
            npw=input(['Desired Nodes per wavelength: (Default ' num2str(npwdef) ') >>']);
        end
        if ~exist('npw','var')||(isempty(npw))
            npw=npwdef;
        end
    end
end



if(zonestrat==1)
    % Nodes per zone
    npzdef=3500;
    if(~usedef) 
        npz=input(['Target Number of Nodes per zone (Default ' num2str(npzdef) ') >>']);
    end
    if ~exist('npz','var')||(isempty(npz))
        npz=npzdef;
    end     
elseif(zonestrat==2)
    %Wavelengths per zone    
    wlperzonedef=0.7;
    if(~usedef) 
        wlperzone=input(['Target Number of wavelengths per zone (Default ' num2str(wlperzonedef) ') >>']);
    end
    if ~exist('wlperzone','var')||(isempty(wlperzone))
        wlperzone=wlperzonedef;
    end 
elseif(zonestrat==3)
    zldef=1.95633e-02;
    if(~usedef) 
        zlength=input(['Zone size (m) (Default ' num2str(zldef) ') >>']);
    end
    if ~exist('wlperzone','var')||(isempty(zlength))
        zlength=zldef;
    end 
end

if(strcmpi(meshtype,'hex'))
    if(meshstrat==3)
        resdef=[0.002 0.002 0.002];
        if(~usedef) 
            targetres=input(['Target Resolution (Default ' num2str(resdef) ') >>']);
        end
        if ~exist('targetres','var')||(isempty(targetres))
            targetres=resdef;
        end

    end
end

% output file stem
%outstmdef='MREv7mesh';
junk=pwd;
dirloc=find((junk=='/')|(junk=='\'));
outstmdef=junk(dirloc(end)+1:end);
%outstmdef='MREv7mesh';

if(strcmpi(meshtype,'hex'))
    if(meshstrat==1)
        outstmdef=[outstmdef '_voxelmesh'];
    elseif(meshstrat==2)
        outstmdef=[outstmdef '_npw' int2str(npw)];
    elseif(meshstrat==3)
        outstmdef=[outstmdef '_res' sprintf('%2.1f-%2.1f-%2.1f',targetres.*1000)];
    end
elseif(strcmpi(meshtype,'tet'))  
    outstmdef=[outstmdef '_tet' num2str(tetmeshres) 'mm']; 
end
    
    
if(~usedef) 
    outstm=input(['output file stem (default is ' outstmdef '):  >>'],'s');
end
if ~exist('outstm','var')||isempty(outstm)
    outstm=outstmdef;
end




if(par)
    matlabpool('open',nlabs);
end

% Displacement Scaling - Many MRE regularization techniques are sensitive
% to the size of the displacements. Either the regularization weights can
% be altered for each case, or the displacements can be scaled so that they
% are always almost the same size.
dispscale = true; % Switch to turn on displacement scaling
dispscalar = 1e-3; % Average displacement amplitude is scaled to this size.
t0=tic;
%% Scale Data if appropriate
A=abs(Ur+1i*Ui);
meanA = mean(A(repmat(mask,[1 1 1 3])==1));
disp(['Mean Displacement Amplitude all directions ' sprintf('%10.3e',meanA)])
if(dispscale)
    disp(['Scaling Displacements to an average size of ' sprintf('%10.3e',dispscalar) 'm'])
    if(meanA==0)
        error('Mean amplitude = 0, cant scale')
    end
    Ur=Ur./meanA.*dispscalar;
    Ui=Ui./meanA.*dispscalar;    
end
%disp(['Check: mean(A(mask)) = ' num2str( mean(A(repmat(mask,[1 1 1 3])==1)) )])

% % Apply a masked selective median filter to clear up outliers
% thresh=0.2; % Thereshold to apply median filter
% disp(['Median filter with threshold ' num2str(thresh)])
% for ii=1:3
%     Ur(:,:,:,ii)=selectivemedianfilter(Ur(:,:,:,ii),mask,thresh);
%     Ui(:,:,:,ii)=selectivemedianfilter(Ui(:,:,:,ii),mask,thresh);
% end



t1=tic;
dim=size(mask);

% Load segmentation file if supplied
if(strcmp('none',regf))
    disp('Using no soft prior segmentation')
    regs=zeros(dim);
    noreg=true;
else
    disp(['Loading segmentation file ' regf])
    load(regf)
    noreg=false;
end

% Check to see if AnatomicalMRI is supplied. Use regs or mask if it is not.
if(~exist('AnatomicalMRI','var'))
    if(noreg)
        disp('Anatomical MRI not supplied, using Mask')
        AnatomicalMRI=mask;
    else
        disp('Anatomical MRI not supplied, using regs')
        AnatomicalMRI=regs;
    end
end
        


if(strcmpi(meshtype,'hex')) %%%% HEXAHEDRAL MESHING
    
    % Apply nans outside mask to use NaN interpolation and not use values
    % outside mask.
    Ur(repmat(mask,[1 1 1 3])==0)=nan;
    Ui(repmat(mask,[1 1 1 3])==0)=nan;
    
    if(meshstrat==1) % Each MR voxel is a FE node
        maskint=mask;
        Urint=Ur;
        clear Ur
        Uiint=Ui;
        clear Ui
        AnatomicalMRI_int=AnatomicalMRI;
        clear AnatomicalMRI
        xin=0:vox(1):vox(1)*(dim(1)-1);
        yin=0:vox(2):vox(2)*(dim(2)-1);
        zin=0:vox(3):vox(3)*(dim(3)-1);
        xout=xin;
        yout=yin;
        zout=zin;   
        regs_int=regs;
        meshres=vox;
    elseif(meshstrat==2)||(meshstrat==3)  %  Process geometry and interpolate to get desired nodes per wavelength if meshstrat = 2.

        if(meshstrat)==2
            targetres(1:3) = wl/npw;
        elseif(meshstrat==3)
           % targetres already defined from inputs;
        end
        disp(['Target Resolution: ' num2str(targetres)])


        % Allow some deviation from the target resolution to fit an integer
        % number of Hex27 elements across the geometry. This will eliminate
        % wasting planes of data when there is an even number of slices which
        % occurs using a voxel based Hex27 meshing scheme.

        Imask=find(mask);
        [mi,mj,mk]=ind2sub(dim,Imask);
        ext(1)=(max(mi)-min(mi))*vox(1);
        ext(2)=(max(mj)-min(mj))*vox(2);
        ext(3)=(max(mk)-min(mk))*vox(3);

        % Number of targetres fitting across ext = ext(ii)/targetres. Want an
        % EVEN number of resolutions across the domain to fit 27 node elements
        % (Even number of resolutions = odd number of nodes) 
        meshres=zeros(1,3);
        for ii=1:3        
            nres=floor(ext(ii)/targetres(ii));
            if(mod(nres,2)==1)
                nres=nres+1;
            end
            meshres(ii) = ext(ii)/nres;
        end
        disp(['Actual New Resolution(m): ' num2str(meshres)])
        disp(['Data Resolution(m): ' num2str(vox)])

        % Warning if mesh is being interpolated to a lower resolution than the
        % data, this is not really a good idea.
        if(meshres(1)>vox(1)||meshres(2)>vox(2)||meshres(3)>vox(3))
            disp([' >> Warning: Mesh is being interpolated to a lower resolution than the data'])        
        end

        % Perform Interpolation
        xin=0:vox(1):vox(1)*(dim(1)-1);
        yin=0:vox(2):vox(2)*(dim(2)-1);
        zin=0:vox(3):vox(3)*(dim(3)-1);


        xout_tmp= (min(mi)-1)*vox(1):meshres(1):(max(mi)-1)*vox(1);
        yout_tmp= (min(mj)-1)*vox(2):meshres(2):(max(mj)-1)*vox(2);
        zout_tmp= (min(mk)-1)*vox(3):meshres(3):(max(mk)-1)*vox(3);

        % Include a buffer around the extents of the mask so that AnatomicalMRI can be
        % interpolated, and the dataset extends past the boundaries of the
        % mask. e.g. for brain data, the buffer means the interpolated AnatomicalMRI
        % shows the skull etc, not just the masked region of the brain.
        bufsiz=[4 4 0]; % Size of buffer in each direction (interpolated MR voxels) MUST BE EVEN
        if(sum(mod(bufsiz,2)>0))
            error('bufsiz must be even to fit hex27 elements efficiently inside mask extents!!')
        end

        xout=zeros(1,2*bufsiz(1)+length(xout_tmp));
        xout(bufsiz(1)+1:bufsiz(1)+length(xout_tmp))=xout_tmp;
        for ii=1:bufsiz(1)
            xout(bufsiz(1)+1-ii)=xout(bufsiz(1)+1)-meshres(1)*ii;
            xout(bufsiz(1)+length(xout_tmp)+ii)=xout(bufsiz(1)+length(xout_tmp))+meshres(1)*ii;
        end

        yout=zeros(1,2*bufsiz(2)+length(yout_tmp));
        yout(bufsiz(2)+1:bufsiz(2)+length(yout_tmp))=yout_tmp;
        for ii=1:bufsiz(2)
            yout(bufsiz(2)+1-ii)=yout(bufsiz(2)+1)-meshres(2)*ii;
            yout(bufsiz(2)+length(yout_tmp)+ii)=yout(bufsiz(2)+length(yout_tmp))+meshres(2)*ii;
        end

        zout=zeros(1,2*bufsiz(3)+length(zout_tmp));
        zout(bufsiz(3)+1:bufsiz(3)+length(zout_tmp))=zout_tmp;
        for ii=1:bufsiz(3)
            zout(bufsiz(3)+1-ii)=zout(bufsiz(3)+1)-meshres(3)*ii;
            zout(bufsiz(3)+length(zout_tmp)+ii)=zout(bufsiz(3)+length(zout_tmp))+meshres(3)*ii;
        end

        if(0==1) % Check locations of xin and xout
            figure
            vjunk=zeros(size(xin));vjunk(min(mi):max(mi))=1;
            plot(xin,vjunk,'r-',xin,zeros(size(xin)),'ro',xout,zeros(size(xout)),'b.')

            xout(1)
            xin(min(mi))

            xout(end)
            xin(max(mi))
            xout
        end

        % Interpolate Ur, Ui and AnatomicalMRI with the function SplineInterp3D_withnans. 
        Urint=zeros([length(xout) length(yout) length(zout) 3]);
        Uiint=zeros([length(xout) length(yout) length(zout) 3]);
        if(par)
            n_int=8;
            needsint=zeros([size(AnatomicalMRI) n_int]);
            needsint(:,:,:,1:3)=Ur;
            needsint(:,:,:,4:6)=Ui;
            needsint(:,:,:,7)=AnatomicalMRI;
            needsint(:,:,:,8)=regs;

            intdata=zeros([length(xout) length(yout) length(zout) n_int]);

            parfor ii=1:n_int
                intdata(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,needsint(:,:,:,ii),xout,yout,zout,2,'no');
            end
            Urint(:,:,:,:)=intdata(:,:,:,1:3);
            Uiint(:,:,:,:)=intdata(:,:,:,4:6);
            AnatomicalMRI_int=intdata(:,:,:,7);
            %regs_int=intdata(:,:,:,8);
            [Yq,Xq,Zq]=meshgrid(yout,xout,zout);  %[Y,X,Z]=meshgrid(4:5,1:3,7:11) 
            regs_int = interp3(yin,xin,zin,regs,Yq,Xq,Zq,'nearest');
            clear intdata needsint Yq Xq Zq
               
        else
            for ii=1:3
                Urint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ur(:,:,:,ii),xout,yout,zout,2,'no');        
            end
            for ii=1:3
                Uiint(:,:,:,ii)=SplineInterp3D_withnans(xin,yin,zin,Ui(:,:,:,ii),xout,yout,zout,2,'no');
            end
            AnatomicalMRI_int=SplineInterp3D_withnans(xin,yin,zin,AnatomicalMRI,xout,yout,zout,2,'no');

            if(noreg)
                regs_int=zeros(size(AnatomicalMRI_int));
            else
                %regs_int=SplineInterp3D_withnans(xin,yin,zin,regs,xout,yout,zout,2,'no');
                [Yq,Xq,Zq]=meshgrid(yout,xout,zout);  %[Y,X,Z]=meshgrid(4:5,1:3,7:11) 
                regs_int = interp3(yin,xin,zin,regs,Yq,Xq,Zq,'nearest');
                clear intdata needsint Yq Xq Zq
            end
        end

    
        clear Ur Ui
        disp('Displacements Interpolated')
        maskint=~isnan(Urint(:,:,:,1));
        clear AnatomicalMRI
        disp('AnatomicalMRI Interpolated')
        if(0==1)
            for ii=1:3
                %montagestack(Ur(:,:,:,ii))
                montagestack(Urint(:,:,:,ii))
            end
            montagestack(AnatomicalMRI_int)
        end    

    end

    if(par)
        matlabpool('close');
    end

    disp('Displacement Processing Complete, Beginning FE Meshing Process')
    tdisp=toc(t1);
    tmesh=tic;

    %Start the meshing process:
    dim=size(maskint);

    % Assign node numbers to each interpolated voxel
    nodnum=1:numel(Urint(:,:,:,1));
    nodnum=reshape(nodnum,size(Urint(:,:,:,1)));

    % Start in bottom corner, march elements along mesh, inserting them if they
    % are inside the mask
    intmp=zeros(prod(floor(dim/2)),27);
    nel=0;
    nodin=false(size(nodnum));
    for ii=1:2:dim(1)-2
        for jj=1:2:dim(2)-2
            for kk=1:2:dim(3)-2
                if(sum(sum(sum(maskint(ii:ii+2,jj:jj+2,kk:kk+2))))==27) % Element is inside mask
                    nel=nel+1;
                    [intmp(nel,:)]=Hex27incidencelist(nodnum(ii:ii+2,jj:jj+2,kk:kk+2)); % Add element to mesh                
                    nodin(ii:ii+2,jj:jj+2,kk:kk+2)=true; % Tag these nodes as included                
                end
            end
        end
    end

    intmp=intmp(1:nel,:);
    tin=toc(tmesh);
    tnodes=tic;
    % Renumber Nodes (exclude unused nodes)
    nod=zeros(prod(dim),3);
    uvwr=zeros(prod(dim),3);
    uvwi=zeros(prod(dim),3);
    AnatomicalMRInod=zeros(prod(dim),1);
    regnod=zeros(prod(dim),1);
    idx=zeros(prod(dim),3);
    old2newnod=zeros(size(nodnum));
    %in=zeros(size(intmp));
    nn=0;
    nbnod=0;

    bnod=zeros(prod(dim),1);

    for ii=1:dim(1)
        for jj=1:dim(2)
            for kk=1:dim(3)
                if(nodin(ii,jj,kk))
                    nn=nn+1;
                    old2newnod(ii,jj,kk)=nn;
                    %in(intmp==nodnum(ii,jj,kk))=nn;
                    nod(nn,:)=[xout(ii) yout(jj) zout(kk)];
                    uvwr(nn,:)=[Urint(ii,jj,kk,1) Urint(ii,jj,kk,2) Urint(ii,jj,kk,3)];
                    uvwi(nn,:)=[Uiint(ii,jj,kk,1) Uiint(ii,jj,kk,2) Uiint(ii,jj,kk,3)];
                    idx(nn,:)=[ii jj kk];
                    AnatomicalMRInod(nn)=AnatomicalMRI_int(ii,jj,kk);
                    regnod(nn)=regs_int(ii,jj,kk);
                end
            end
        end    
    end
    in=old2newnod(intmp); % 330 times faster than in(intmp==nodnum(ii,jj,kk))=nn; inserted in nested loop above. Checked output files are the same with 'diff' command.
    clear intmp
    tnodrenum=toc(tnodes);

    % Resize arrays to correct size.
    nod=nod(1:nn,:);
    uvwr=uvwr(1:nn,:);
    uvwi=uvwi(1:nn,:);
    AnatomicalMRInod=AnatomicalMRInod(1:nn);
    regnod=regnod(1:nn);
    idx=idx(1:nn,:);
    
    tbni=tic;
    disp('Building Element Connetivity')
    % Build element connectivity number to determine boundary nodes
    elcon=zeros(nn,1);
    for ii=1:nel
        for jj=1:27
            elcon(in(ii,jj))=elcon(in(ii,jj))+1;
        end
    end
    disp('Finding Boundary nodes')
    disp(['nn = ' int2str(nn)])
    % Build Boundary node array based on nodal connectivity % THis step accounts for ~99% of the meshing time 
    bnum=1;
    %New bnode finding - 2800x faster than old version - bnodes come out
    %in differet order but get all the same nodes.
    bnodchecked=false([nn 1]);
    for ii=1:nel
        for jj=1:27
           if(~bnodchecked(in(ii,jj)))
                if jj<=4 || (jj>=10 && jj<=13) % corner nodes =1,2,3,4,10,11,12,13
                    if elcon(in(ii,jj))<8 %corner node is boundary node
                        bnod(bnum)=in(ii,jj);
                        bnum=bnum+1;
                    end
                elseif (jj>=5 && jj<=8) || (jj>=14 && jj<=17) || (jj>=19 && jj<=22) % mid-edge nodes = 5,6,7,8,14,15,16,17,19,20,21,22
                    if elcon(in(ii,jj))<4 %mid-edge node is boundary node
                        bnod(bnum)=in(ii,jj);
                        bnum=bnum+1;
                    end
                elseif  jj==9 || jj==18 || (jj>=23 && jj<=26) %mid-face nodes are 9,18,23,24,25,26
                    if elcon(in(ii,jj))<2 %mid-face node is boundary node
                        bnod(bnum)=in(ii,jj);
                        bnum=bnum+1;
                    end
                end
                bnodchecked(in(ii,jj))=true;
           end               
        end
    end


    nbnod=bnum-1;
    bnod=bnod(1:nbnod);
    tbnod=toc(tbni);

    tfemesh=toc(tmesh); % Time for full meshing process

    disp(['Hexahedral Meshing Complete, ' int2str(nn) ' Nodes and ' int2str(nel) ' elements'])

elseif(strcmpi(meshtype,'tet'))  %%%% TETRAHEDRAL MESHING
    tetmeshres=tetmeshres*0.85; % Adjust mesh resolution as bgmesh creates 
                                % elements that are slightly bigger than the specified resolution
    
    %% Generate input file for tet meshing (MREmesh.m)
    % Mesh locations
    xin=0:vox(1):vox(1)*(dim(1)-1);
    yin=0:vox(2):vox(2)*(dim(2)-1);
    zin=0:vox(3):vox(3)*(dim(3)-1);
    
    spmeshfile=false;
    if(spmeshfile)
        fid=fopen('Tetmeshinginput','w');
        fprintf(fid,['3 \n']);
        fprintf(fid,[outstm '\n']);
        fprintf(fid,[num2str(tetmeshres/1000) ' \n']);
        fprintf(fid,['YOUREMAILADDRESS@dartmouth.edu \n']);
        fprintf(fid,['0.85d0 \n']);
        fprintf(fid,['100 \n']);
        fprintf(fid,['1d-4 \n']);
        fprintf(fid,['16 \n']);
        fclose(fid);
        disp('TET MESH SPECIFIED: ')
        disp('Run MREmesh(''Tetmeshinginput'')')
        disp(['Move all mesh files to tet/' outstm ' directory'])
        disp('Command to put column of ones in node file :: ')
        disp('awk ''{print $1,$2,$3,$4,1}'' meshstem.nod > meshstem.hom.nod')
        disp('awk ''{print $1,$2,$3,$4,1}'' meshstem.nod > meshstem.hom.nod')
        disp('Create pressure BC file using')
        disp('addporopressureBConds_v8.m')
        disp('Update line 11 of poro runfile to correct filename')
        disp('It would be good to have some better meshing software here')
    end
    
    %% Tetrahedral meshing taken from Ligins code
    stackfn=logical(mask);
    save stackfn.mat stackfn;

    slicethickness=vox(3)*1000; 

    stackInfo=struct('Rows',            mridim(1),...
                      'Columns',        mridim(2),...
                      'PixelSpacing',   [vox(1) vox(2)]*1000,...
                      'SliceThickness', slicethickness);
                      %'SpacingBetweenSlices',SpacingBetweenSlices);

    save stackInfo.mat stackInfo;

    %% Create Surface mesh  
    
    isoval=0.99; % The level that the isosurface is set at. If it was set 
                 % at 1, the top and bottom slices were on the very edge 
                 % of the data which caused interpolation issues when 
                 % interpolating back onto the data.     
    tic
    [eb, pb]=CreateSurfaceFrom3DImage('stackfn.mat','stackInfo.mat',tetmeshres, isoval);
    
    if(strcmp(smoothmesh,'y'))
        % Smooth mesh surface (horribly inefficent)
        disp('Beginning Horribly inefficent mesh smoothing')
        sigma=tetmeshres*3;
        cutrad=3*sigma;
        pbsmooth=zeros(size(pb));
        wsum=zeros(size(pb,1),1);
        for ii=1:length(pb)
            for jj=1:length(pb)
                dij=sqrt((pb(ii,1)-pb(jj,1))^2+(pb(ii,2)-pb(jj,2))^2+(pb(ii,3)-pb(jj,3))^2);
                if(dij<cutrad)
                    gaussfac=exp(-0.5*(dij/sigma)^2);
                    pbsmooth(ii,1)=pbsmooth(ii,1)+gaussfac*pb(jj,1);
                    pbsmooth(ii,2)=pbsmooth(ii,2)+gaussfac*pb(jj,2);
                    pbsmooth(ii,3)=pbsmooth(ii,3)+gaussfac*pb(jj,3);
                    wsum(ii)=wsum(ii)+gaussfac;
                end
            end
        end
        pb(:,1)=pbsmooth(:,1)./wsum;
        pb(:,2)=pbsmooth(:,2)./wsum;
        pb(:,3)=pbsmooth(:,3)./wsum;
        
        disp('Horrifically inefficent smoothing complete')
    end
    
    save surfacemesh.mat eb pb

    % Plot mesh to check
    fig=figure;
    fig.Position=[200 700 1100 400];
    subplot (1,2,1)
    plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
    AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);
    axis equal; axis(AXIS);
    grid on;
    title('Surface Mesh Nodes'); hold on;    
        
    %% Start of Internal Mesh
    if tetcustom % Customize mesh criteria
        facet_angle=input('Facet angle (default=25)>>: '); 
        facet_size=input(['Facet size (default= ' num2str(tetmeshres) ')>>:']);
        facet_distance=input('Facet distance (default=0.6)>>:');
        cell_radius_edge=input('Cell radius edge (default=2)>>:');
        cell_size=input(['Cell size (default=' num2str(tetmeshres) ')>>:']);
        tumor_cen=input('Tumor Center (default=[100.0000 100.000 20.00])>>:');
        surf_cen=input('Surface Center (default=[100.0000 100.0000 20.000];)>>:');
        R1=input('R1 (default=30)>>:');
        S1=input(['S1 (default=' num2str(tetmeshres) ')>>:']);
        R2=input('R2 (default=5)>>:');
        S2=input(['S2 (default=' num2str(tetmeshres) ')>>:']);
        ref_ratio=input('Reference Ratio (default=1)>>:');
    else
        facet_angle=[25]; 
        facet_size=[tetmeshres];
        facet_distance=[0.6];
        cell_radius_edge=[2];
        cell_size=[tetmeshres];
        tumor_cen=[100.0000 100.000 20.00]; % We dont define a 'tumor', so keep these cen values equal
        surf_cen=[100.0000 100.0000 20.000];
        R1=[30];  % Not entirerly sure what these do
        S1=[tetmeshres];
        R2=[5];
        S2=[tetmeshres];
        ref_ratio=[1];        
    end
    criteria=struct('facet_angle',facet_angle,'facet_size',facet_size, ...
                    'facet_distance',facet_distance, 'cell_radius_edge',cell_radius_edge, ...
                    'cell_size', cell_size, 'tumor_cen',tumor_cen, 'surf_cen',surf_cen, ...
                    'R1',R1, 'S1',S1, 'R2',R2, 'S2',S2', 'ref_ratio',ref_ratio);

    disp(['Surface Mesh created: ' int2str(round(toc)) ' seconds']);
    %% Main meshing step
    tic;
    if(exist('~/outmesh.mesh','file'))
      delete('~/outmesh.mesh');
    end
    [in,nodtemp]=MakeTetraMesh(eb,pb,criteria);
    disp(['Tetrahedral Meshing Complete ' int2str(round(toc)) ' seconds'])
    

    %Switched the nodes to match the hexahedral meshing!-LS 01/19/16
    nod(:,2)=nodtemp(:,1)*1e-3; 
    nod(:,1)=nodtemp(:,2)*1e-3; 
    if min(nodtemp(:,3))<0
        nod(:,3)=nodtemp(:,3)*1e-3+min(nodtemp(:,3))*1e-3;
        warning('Nodes with negative depth appeared. Should be no problem.....')
    else
        nod(:,3)=nodtemp(:,3)*1e-3;
    end
    
    %% Generating Boundary Nodes and Elements
    [bel]=getBdyFromMesh(in(:,1:4)); 
    bnod=unique(bel); 
    
    % Plot all nodes and boundary nodes
    subplot(1,2,2);
    fig.Position=[200 700 1100 1000];
    plot3(nod(~bnod,1),nod(~bnod,2),nod(~bnod,3),'b.')
    hold on;
    plot3(nod(bnod,1),nod(bnod,2),nod(bnod,3),'r.');
    title('All Nodes');
    xlabel('X'); ylabel('Y'); zlabel('Z')
    grid on; AXIS = ([0 max(nod(:,1)) 0 max(nod(:,2)) 0 max(nod(:,3))]);
    axis equal; axis(AXIS);
    legend('Interior Nodes','Boundary Nodes','Location','Best')

    %% Sanity checks
    if isempty(in); 
        warning('There are no elements! Remesh using different criteria.');
    end        
    if isempty(bel)
        warning('There does not seem to be any boundary elements')
    end
    if isempty(bnod);
        warning('There does not seem to be any boundary nodes')
    end

    %% Interpolate displacements from grid to the FE mesh    
    disp('Interpolating Motion Data')
    [uvwr,uvwi]=tet_interp_disp(nod,Ur,Ui,vox,mask,AnatomicalMRI);
    
    %% Interpolate regions to the FE mesh
    regnod=tet_nearest_interp(nod,regs,vox);
    
    %% Interpolate AnatomicalMRI to the FE mesh
    AnatomicalMRInod=tet_nearest_interp(nod,AnatomicalMRI,vox);
    
    
    disp('Mesh Created!');

    disp(sprintf('Number of Boundary Nodes: %d \r',size(bnod,1))); 
    disp(sprintf('Number of Nodes: %d \r',size(nod,1)));
    disp(sprintf('Number of Internal Nodes: %d \r',size(nod,1)-size(bnod,1))) 
    disp(['Number of elements ' int2str(length(in))])        
end

% Find bottom and top surface nodes for poro pressure BCs
rngz=[min(nod(:,3)) max(nod(:,3))];
btmbnods=(nod(bnod,3)<(rngz(1)+0.02*diff(rngz)));
topbnods=(nod(bnod,3)>(rngz(2)-0.02*diff(rngz)));   


% Perform coordinate transformation
uvwr=uvwr*Motion2Image';
uvwi=uvwi*Motion2Image';  % Nnx3*3x3=Nn*3, Transpose IS required because we store uvw as COLUMN vectors, and transformation matrices are defined for ROW vectors [x' y' z']=T[x y z]

% Reverse z direction if this is a left handed coordinate system
if(~RHcoord)
    disp('RHcoord indicates Coordinate system not Right handed. Reversing z direction')
    nod(:,3)=-nod(:,3);
    uvwr(:,3)=-uvwr(:,3);
    uvwi(:,3)=-uvwi(:,3);
end
    

%% Calcualte Zone sizes
touti=tic;
% Zonestrat=1: target number of nodes per zone
% Zonestrat=2: target number of wavelengths per zone
% Zonestrat=3: Define subzone size
znovlp=[0.15 0.15 0.15]; % Zone overlap
znovlp1=1+znovlp;

if(zonestrat==1)
    
    % -> Aim for approximatly cubic zones, with a specified number of nodes per zone
    % npz = Target number of nodes per zone (Note that this strategy usually ends up with zones ~70% of target size).
    disp(['zonestrat=1, aiming for ' int2str(npz) ' nodes per zone'])
    zfac=rangedim./min(rangedim); % Ratio of mesh dimensions
   
    Nz=nn*znovlp1(1)*znovlp1(2)*znovlp1(3)/npz;
    F=(Nz/(zfac(1)*zfac(2)*zfac(3)))^(1/3);
    znedge=zfac.*F;
   
    % v7.05+ = direct specification of zone sizes, [L1 L2 L3]
    % Aim for cubic zones, NPZ =(L*(1+2*ovlp(1)))/meshres(1))*(L*(1+2*ovlp(2)))/meshres(2))*(L*(1+2*ovlp(3)))/meshres(3))
    % L^3=NPZ*mesres(1)*meshres(2)*meshres(3)/(1+2*ovlp(1))*(1+2*ovlp(2))*1+2*ovlp(3)))
    znedgelength(1:3)=( npz*meshres(1)*meshres(2)*meshres(3) / ((1+2*znovlp(1))*(1+2*znovlp(2))*(1+2*znovlp(3))) )^(1/3);
    disp(['Zone edge length ' num2str(znedgelength)])
    disp(['Estimated nodes per zone = ' int2str(prod(znedgelength.*(1+2*znovlp)./meshres))])
elseif(zonestrat==2)
    disp(['zonestrat=2, aiming for ' int2str(wlperzone) ' wavelengths per zone'])    
     
    znesz=wlperzone*wl;

    % Actual SZ length = (1+2*znovlp(i))*L(i)
    for ii=1:3
        znedgelength(ii)=znesz/(1+2*znovlp(ii));
    end
    disp(['Zone edge length ' num2str(znedgelength)])
    disp(['Estimated nodes per zone = ' num2str(prod(znesz./meshres)) ])
elseif(zonestrat==3)
    
    for ii=1:3
        znedgelength(ii)=zlength;
    end
    disp(['Zone edge length ' num2str(znedgelength)])
end



%% Output mesh
% Create Directories if they dont exist
D=dir; % Check to see if hex exists
direc=pwd;
if strcmpi(meshtype,'hex')
    hexflg=0;
    for ii=1:length(D)
        if length(D(ii).name)==3
            if (strcmp(D(ii).name,'hex')&&(D(ii).isdir==1)) 
                hexflg=1;
            end
        end
    end
    if (hexflg==0) % Make Hex
        disp('Creating hex directory')
        mkdir('hex')
    end
    junk=['hex/' outstm]; % Make subdirectiories
    mkdir('hex',outstm);
    mkdir(junk,'inv');
    inpath=[direc '/hex/' outstm '/'];        
elseif(strcmpi(meshtype,'tet'))    
    tetflg=0;
    for ii=1:length(D)
        if length(D(ii).name)==3
            if (strcmp(D(ii).name,'tet')&&(D(ii).isdir==1)) 
                tetflg=1;
            end
        end
    end
    if (tetflg==0) % Make tet
        disp('Creating tet directory')
        mkdir('tet')
    end
    junk=['tet/' outstm]; % Make subdirectiories
    mkdir('tet',outstm);
    mkdir(junk,'inv');
    inpath=[direc '/tet/' outstm '/'];
    % Move Meshing files into the output mesh directory
    movefile('surfacemesh.mat',inpath);
    movefile('stackInfo.mat',inpath);
    movefile('stackfn.mat',inpath);
end

outpath=['inv/'];

%mask file
maskoutf=[outstm '.mask.mat'];
save([inpath maskoutf],'mask');

%nod files
ind=(1:size(nod,1))';
nodhmgf=[outstm '.hom.nod'];
nodregoutf=[outstm '.reg.nod'];

fid=fopen([inpath nodhmgf],'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e  1 \n',[ind nod]');
fclose(fid);

fid=fopen([inpath nodregoutf],'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %3i \n',[ind nod regnod]');
fclose(fid);

elmind=(1:size(in,1))';
if(strcmpi(meshtype,'hex'))  % Output hexahedral mesh files
   
    % Element file  (extra column of ones to indicate
    % homogeneous elementally defined properties
    elmhmgf=[outstm '.hom.elm'];
    fid=fopen([inpath elmhmgf],'w');
    fprintf(fid,'%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i  1 \n',[elmind in]');
    fclose(fid);


    %index file -> Note this is interpolated index. Need to back-interpolate
    %using xin and xout to get values at original MR voxels.
    idxf=[outstm '.idx'];
    fid=fopen([inpath idxf],'w');
    fprintf(fid,'%7i %6i %6i %6i\n',[ind idx]');
    fclose(fid);
    
elseif(strcmpi(meshtype,'tet'))
    % Element file  (extra column of ones to indicate
    % homogeneous elementally defined properties
    elmhmgf=[outstm '.hom.elm'];
    fid=fopen([inpath elmhmgf],'w');
    fprintf(fid,'%8i %8i %8i %8i %8i %8i \n',[elmind in]');
    fclose(fid);
    
    
    % Create some pressure BC files
    pbcf1=[outstm 'typ1.pbcs'];
    fid=fopen([inpath pbcf1],'w');
    fprintf(fid,'%8i %8i %8i %15.6e %15.6e \n',[(1:length(bnod))' bnod ones(length(bnod),1) zeros(length(bnod),1) zeros(length(bnod),1)]');
    fclose(fid);    
    
    typ=2*ones(size(bnod));
    typ(btmbnods)=1;
    typ(topbnods)=1;
    
    pbcf2=[outstm 'typ2sides.pbcs'];
    fid=fopen([inpath pbcf2],'w');
    fprintf(fid,'%8i %8i %8i %15.6e %15.6e \n',[(1:length(bnod))' bnod typ zeros(length(bnod),1) zeros(length(bnod),1)]');
    fclose(fid);    
    
    typ=ones(size(bnod));
    typ(btmbnods)=2;
        
    pbcf3=[outstm 'typ2btm.pbcs'];
    fid=fopen([inpath pbcf3],'w');
    fprintf(fid,'%8i %8i %8i %15.6e %15.6e \n',[(1:length(bnod))' bnod typ zeros(length(bnod),1) zeros(length(bnod),1)]');
    fclose(fid);    
    
    % Boundary element file
    belf=[outstm '.bel'];
    fid=fopen([inpath belf],'w');
    fprintf(fid,'%8i %8i %8i %8i \n',[(1:length(bel))' bel]');
    fclose(fid); 
        
end

    
%disp file
dspoutf=[outstm '.dsp'];
fid=fopen([inpath dspoutf],'w');
fprintf(fid,'%7i %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e \n',[ind uvwr(:,1) uvwi(:,1) uvwr(:,2) uvwi(:,2) uvwr(:,3) uvwi(:,3)]');
fclose(fid);

%anatomical.mtrin file
magoutf=[outstm 'anatomical.mtrin'];
fid=fopen([inpath magoutf],'w');
fprintf(fid,'%7i %15.8e\n',[ind AnatomicalMRInod]');
fclose(fid);


%bnod file
bcoutf=[outstm '.bnd'];
fid=fopen([inpath bcoutf],'w');
fprintf(fid,'%7i  %7i\n',[(1:length(bnod))' bnod]');
fclose(fid);

if(strcmpi(meshtype,'hex'))
    %save data required to get back to original MR voxels (including mask filename)
    if isempty(msk)
        msk=maskoutf;
    end
    save([inpath outstm '.InterpLocations.mat'],'maskint','nodin','xout','yout','zout','xin','yin','zin','msk') 

    %save interpolated motion data for possible later use
    if(meshstrat~=1)
        save([inpath outstm '.InterpData.mat'],'maskint','Urint','Uiint','AnatomicalMRI_int','Motion2Image')
    end
end

%save inputs to meshing code to link with mesh
save([inpath outstm '.meshinput.mat'],'zonestrat','dispscale','dispscalar','meshtype')

if(zonestrat==2)
    save([inpath outstm '.meshinput.mat'],'wlperzone','-append')
end

if(zonestrat==2)
    save([inpath outstm '.meshinput.mat'],'wlperzone','-append')
end

if(strcmpi(meshtype,'tet'))
    save([inpath outstm '.meshinput.mat'],'tetmeshres','-append')
    save([inpath outstm '.meshinput.mat'],'criteria','-append')
end

if(strcmpi(meshtype,'hex'))
    save([inpath outstm '.meshinput.mat'],'meshstrat','-append')   
    if((meshstrat==2)||(zonestrat==2))
        save([inpath outstm '.meshinput.mat'],'muest','rhoest','-append')
    end
    if(meshstrat==2)||(meshstrat==3)
        save([inpath outstm '.meshinput.mat'],'targetres','meshres','-append')
    end
    if(meshstrat==2)
        save([inpath outstm '.meshinput.mat'],'npw','-append')
    end
end

% Output imginfo to allow interpolation back to the original MR voxels
fid=fopen([inpath outstm '.imginfo'],'w');
fprintf(fid,'Size of MRI stack \n');
fprintf(fid,'%7i %7i %7i \n',mridim);
fprintf(fid,'Voxel dimensions (mm) \n');
fprintf(fid,'%15.6e %15.6e %15.6e \n',voxsize_mm);
fprintf(fid,'Right hand coordinate system indicator \n');
fprintf(fid,'%i',RHcoord);
fclose(fid);


% region file for soft prior -> Output the region file at the original data
% resolution so that the regions can be interpolated to the material
% property meshes inside the recon code.
regoutf=[outstm '.regstack'];

writeregstackfile(regs,mridim,xin,yin,zin,[inpath regoutf])

% Write .sfweight file to supply variable out-of-region smoothing weight
sprsf_file=[inpath outstm '.regstack.sfweight'];
fid=fopen(sprsf_file,'w');
fprintf(fid,'%5.2f',0.5*ones([1 max(regs(:))]));
fclose(fid);  


% Load DTI if supplied
% From Curtis: 1st index is LR, 2nd index is AP, 3rd is SI. 

if(strcmp('none',DTIf))
    disp('no DTI file supplied')
    noDTI=true;
else
    disp(['Loading DTI file ' DTIf])
    load(DTIf)
    noDTI=false;
end

if ~noDTI
    DTIstackfx=[inpath outstm '.DTIstackx'];
    DTIstackfy=[inpath outstm '.DTIstacky'];
    DTIstackfz=[inpath outstm '.DTIstackz'];
    
    
    simdata=true;
    
    if(simdata)
        disp('using simulated DTI direction assumption that V1(:,:,:,1) = (:,1,1), V1(:,:,:,2) = (1,:,1), V1(:,:,:,3) = (1,1,:)')
        writeregstackfile(V1(:,:,:,1),size(V1(:,:,:,1)),xin,yin,zin,DTIstackfx); % First index is LR, which is (1,:,1), which is our Y
        writeregstackfile(V1(:,:,:,2),size(V1(:,:,:,2)),xin,yin,zin,DTIstackfy); % 2nd index is AP, which is (:,1,1), which is our X
        writeregstackfile(V1(:,:,:,3),size(V1(:,:,:,3)),xin,yin,zin,DTIstackfz); % 3rd index is SI, which is (1,1,:), which is our Z
        
    else
        disp('using MRI DTI direction assumption that V1(:,:,:,1) = (1,:,1), V1(:,:,:,2) = (:,1,1), V1(:,:,:,3) = (1,1,:)')

        writeregstackfile(V1(:,:,:,1),size(V1(:,:,:,1)),xin,yin,zin,DTIstackfy); % First index is LR, which is (1,:,1), which is our Y
        writeregstackfile(V1(:,:,:,2),size(V1(:,:,:,2)),xin,yin,zin,DTIstackfx); % 2nd index is AP, which is (:,1,1), which is our X
        writeregstackfile(V1(:,:,:,3),size(V1(:,:,:,3)),xin,yin,zin,DTIstackfz); % 3rd index is SI, which is (1,1,:), which is our Z
    end
end
%% Output Runfiles and submit files




if(strcmpi(meshtype,'hex')) % Output Hexahedral runfiles
    output_runfiles_hex_v9_liang(outstm,vox,muest,rhoest,DR,inpath,nodhmgf,elmhmgf,bcoutf,regoutf,freqHz,dspoutf,outpath,znedgelength,znovlp,noreg,noDTI,SubjectName);
elseif(strcmpi(meshtype,'tet'))
    output_runfiles_tet_v9_liang(outstm,vox,muest,rhoest,DR,inpath,nodhmgf,elmhmgf,bcoutf,regoutf,pbcf2,freqHz,dspoutf,outpath,znedgelength,znovlp,noreg,noDTI,SubjectName);
end

% %% Move all files to where they should be.
% for ii=1:length(mvfiles)
% movefile(mvfiles(ii).name,inpath); 
% end
tout=toc(touti);
ttotal=toc(t0);
%% Disply time for each part:

if(strcmp(meshtype,'hex'))
disp(['Time for Displacment processing: ' sprintf('%6.2f',tdisp) ' seconds'])
disp(['Time for FE meshing Process: ' sprintf('%6.2f',tfemesh) ' seconds'])
disp(['  Time to build incidence list: ' sprintf('%6.2f',tin) ' seconds'])
disp(['  Time to renumber nodes: ' sprintf('%6.2f',tnodrenum) ' seconds'])
disp(['  Time to find boundary nodes: ' sprintf('%6.2f',tbnod) ' seconds'])
disp(['Time to output files: ' sprintf('%6.2f',tout) ' seconds'])
disp(' ')
disp(['Total processing time: ' sprintf('%6.2f',ttotal) ' seconds'])
end
end

function [in]=Hex27incidencelist(nodnum)
%Hex27incidencelist:outputs the row of the incidence list from a Hex27
%element built of of a 3x3x3 cube of node numbers, nodnum.

% Nodal coords copied directly from matts masters thesis.
% X(1,:) = [-1,-1,-1 ]; X(10,:) = [-1,-1, 1 ]; X(19,:) = [-1,-1, 0 ];
% X(2,:) = [ 1,-1,-1 ]; X(11,:) = [ 1,-1, 1 ]; X(20,:) = [ 1,-1, 0 ];
% X(3,:) = [ 1, 1,-1 ]; X(12,:) = [ 1, 1, 1 ]; X(21,:) = [ 1, 1, 0 ];
% X(4,:) = [-1, 1,-1 ]; X(13,:) = [-1, 1, 1 ]; X(22,:) = [-1, 1, 0 ];
% X(5,:) = [ 0,-1,-1 ]; X(14,:) = [ 0,-1, 1 ]; X(23,:) = [ 0,-1, 0 ];
% X(6,:) = [ 1, 0,-1 ]; X(15,:) = [ 1, 0, 1 ]; X(24,:) = [ 1, 0, 0 ];
% X(7,:) = [ 0, 1,-1 ]; X(16,:) = [ 0, 1, 1 ]; X(25,:) = [ 0, 1, 0 ];
% X(8,:) = [-1, 0,-1 ]; X(17,:) = [-1, 0, 1 ]; X(26,:) = [-1, 0, 0 ];
% X(9,:) = [ 0, 0,-1 ]; X(18,:) = [ 0, 0, 1 ]; X(27,:) = [ 0, 0, 0 ];
% Plot the nodal coords to check
% for ii=1:27
%     plot3(X(ii,1),X(ii,2),X(ii,3),'r.','markersize',18)
%     title(['node ' int2str(ii) ' added'])
%     grid on
%     hold on
%     pause
% end
% Added 2 to X in matlab to get indices for nodnum 

in=zeros(1,27);

in(1)=nodnum(1,1,1);
in(2)=nodnum(3,1,1);
in(3)=nodnum(3,3,1);
in(4)=nodnum(1,3,1);
in(5)=nodnum(2,1,1);
in(6)=nodnum(3,2,1);
in(7)=nodnum(2,3,1);
in(8)=nodnum(1,2,1);
in(9)=nodnum(2,2,1);
in(10)=nodnum(1,1,3);
in(11)=nodnum(3,1,3);
in(12)=nodnum(3,3,3);
in(13)=nodnum(1,3,3);
in(14)=nodnum(2,1,3);
in(15)=nodnum(3,2,3);
in(16)=nodnum(2,3,3);
in(17)=nodnum(1,2,3);
in(18)=nodnum(2,2,3);
in(19)=nodnum(1,1,2);
in(20)=nodnum(3,1,2);
in(21)=nodnum(3,3,2);
in(22)=nodnum(1,3,2);
in(23)=nodnum(2,1,2);
in(24)=nodnum(3,2,2);
in(25)=nodnum(2,3,2);
in(26)=nodnum(1,2,2);
in(27)=nodnum(2,2,2);
end

%% Thresholded Median filter
function [stackout]=selectivemedianfilter(stackin,mask,thresh)

s=size(stackin);
stackout=stackin;
fsz=1;

for ii=1:s(1)
    %disp(['ii = ' int2str(ii)])
    for jj=1:s(2)
        for kk=1:s(3)
            if(mask(ii,jj,kk)==1)
                nmask=mask( max(ii-fsz,1):min(ii+fsz,s(1)),max(jj-fsz,1):min(jj+fsz,s(2)),max(kk-fsz,1):min(kk+fsz,s(3)) )==1;
                nstack=stackin( max(ii-fsz,1):min(ii+fsz,s(1)),max(jj-fsz,1):min(jj+fsz,s(2)),max(kk-fsz,1):min(kk+fsz,s(3)) );
                medv=median(nstack(nmask));
                if( abs((medv-stackin(ii,jj,kk))/medv)>thresh )
                    stackout(ii,jj,kk)=medv;
                end
            end
        end
    end
end

end

%% Clean hanging mask voxels (i.e. ones that stic out only 1 voxel thick in any direction)
function maskout=cleanmask(mask)

    s=size(mask);
    maskout=zeros(size(mask));
    
    for ii=1:s(1)
        for jj=1:s(2)
            for kk=1:s(3)
                if(mask(ii,jj,kk))
                    d1=sum(mask(max(1,ii-1):min(s(1),ii+1),jj,kk));
                    d2=sum(mask(ii,max(1,jj-1):min(s(2),jj+1),kk));
                    d3=sum(mask(ii,jj,max(1,kk-1):min(s(3),kk+1)));
                    if((d1>1)&&(d2>1)&&(d3>1))
                        maskout(ii,jj,kk)=1;
                    else
                        maskout(ii,jj,kk)=0;
                    end
                end
            end
        end
    end
    
    %montagestack(mask+maskout,[],'clip');colormap(gca,jet)
    %title('Removed Hanging voxels=green')
    

end

function writeregstackfile(stack,stackdim,x1,x2,x3,regoutf)

% Format:
% StackSize
% x coords
% y coords
% z coords
% slice1
% slice2
% ..
% Slice_end
fid=fopen(regoutf,'w');
fprintf(fid,'Array Dimensions \n'); 
fprintf(fid,'%7i %7i %7i \n',stackdim);
fprintf(fid,'x coordinates \n');
fclose(fid);
dlmwrite(regoutf,x1,'-append','delimiter',' '); 
fid=fopen(regoutf,'a'); fprintf(fid,'y coordniates \n'); fclose(fid);
dlmwrite(regoutf,x2,'-append','delimiter',' ');
fid=fopen(regoutf,'a'); fprintf(fid,'z coordniates \n'); fclose(fid);
dlmwrite(regoutf,x3,'-append','delimiter',' ');
for ii=1:stackdim(3)
    dlmwrite(regoutf,stack(:,:,ii),'-append','delimiter',' ');
end


end
            

  

