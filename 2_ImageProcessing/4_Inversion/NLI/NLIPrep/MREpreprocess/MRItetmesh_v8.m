function MRItetmesh_v8(default)

%% This is a tetrahedral meshing program 
% Adapted from bgmesh (Brain Group Meshing with tetrahedral elements)
% Converted from bgmesh.m to MRItetmesh_v8 to be incorporated into
% MREMotion2.m 
% Author: Ligin Solamen
% PhD Candidate
% July 2016
% Thayer School of Engineering

% Lines changed 201-206 for testing boundary nodes 

%% Parallel ToolBox and default mode
par=false;          % If you have the matlab parallel processing toolbox set this to true, otherwise false.
nlabs=7;            % Number of labs for matlab to use (max 7)
%[usedef]=ParDefault(nlabs,par,nargin,default);
disp('MRE-Zone v7.3 Data Converter')
if(par)
    disp(['Using Parallelized version with ' int2str(nlabs) ' labs'])
else
    disp('Using non-Parallelized version - modify value of ''par'' to change')
end

if(nargin==0)       % Default value is to prompt for inputs.
    default='no';
end
usedef=strcmp(default,'default'); % If no => false;else=>true
if(usedef)
    disp('Using Default values, no input prompts')   
end

[default]=SetDefault_v8;

[A,vox,mridim,P,freqHz,msk,mask,MREMotionInfo]=DataExtraction_v8(usedef,default);

[strat]=Strategies_v8(usedef,default);

%% Get Inputs

stackfn=logical(mask);
save stackfn.mat stackfn;

slicethickness=(MREMotionInfo(1).slicegap+MREMotionInfo(1).slicethickness); 

stackInfo=struct('Rows',            MREMotionInfo(1).nX,...
                  'Columns',        MREMotionInfo(1).nY,...
                  'PixelSpacing',   [MREMotionInfo(1).pixelspacing_x MREMotionInfo(1).pixelspacing_y],...
                  'SliceThickness', slicethickness);
                  %'SpacingBetweenSlices',SpacingBetweenSlices);

save stackInfo.mat stackInfo;

%% Create Nodes and Elements

load stackfn.mat stackfn
montagestack(double(stackfn)); drawnow; 
load stackInfo.mat stackInfo

[eb, pb]=CreateSurfaceFrom3DImage('stackfn.mat','stackInfo.mat');

save surfacemesh.mat eb pb


fig=figure;
fig.Position=[200 700 1100 400];
subplot (1,2,1)
%plot3(pb(:,1), pb(:,2), pb(:,3), '.'); 
plot3(pb(:,1)*1e-3, pb(:,2)*1e-3, pb(:,3)*1e-3, '.'); 
AXIS = ([0 max(pb(:,1))*1e-3 0 max(pb(:,2))*1e-3 0 max(pb(:,3))*1e-3]);
%AXIS = ([0 max(pb(:,1)) 0 max(pb(:,2)) 0 max(pb(:,3))]);
axis equal; axis(AXIS);
%'Color',[0.5 0 1]
grid on;
title('Surface Mesh Nodes'); hold on;
fname=input('Name of input file: (Default: input) >>');
if isempty(fname)==1
    fname='input';
end

fid=fopen(fname,'rt');
    

todo = fgetl(fid);
disp(['Selection --> ',todo]);
todo = str2num(todo);
fprintf('\nEnter the series name: \n');
%read series namee
outstmtemp = fgetl(fid);
%disp(['Series Name --> ',outstm]);

outstm=input(strcat('Output Stem (Default :',outstmtemp,')>>'),'s');
if isempty(outstm)==1
    outstm=outstmtemp;
end





%% Start of Internal Mesh
temp=input('Do you want to customize your tetrahedral meshing criteria? 1=yes, (Default: no)>>:');
if temp==1
    facet_angle=input('Facet angle (default=22)>>: '); 
    facet_size=input('Facet size (default=3)>>:');
    facet_distance=input('Facet distance (default=0.8)>>:');
    cell_radius_edge=input('Cell radius edge (default=3)>>:');
    cell_size=input('Cell size (default=3)>>:');
    tumor_cen=input('Tumor Center (default=[100.0000 100.000 20.00])>>:');
    surf_cen=input('Surface Center (default=[100.0000 100.0000 20.000];)>>:');
    R1=input('R1 (default=30)>>:');
    S1=input('S1 (default=3)>>:');
    R2=input('R2 (default=5)>>:');
    S2=input('S2 (default=3)>>:');
    ref_ratio=input('Reference Ratio (default=1)>>:');
    criteria=struct('facet_angle',facet_angle,'facet_size',facet_size, ...
                    'facet_distance',facet_distance, 'cell_radius_edge',cell_radius_edge, ...
                    'cell_size', cell_size, 'tumor_cen',tumor_cen, 'surf_cen',surf_cen, ...
                    'R1',R1, 'S1',S1, 'R2',R2, 'S2',S2', 'ref_ratio',ref_ratio);
                
    
    if isempty(criteria.facet_angle)
        criteria.facet_angle=22;
    end
    
    if isempty(criteria.facet_size)
        criteria.facet_size=3;
    end
    
    if isempty(criteria.facet_distance)
        criteria.facet_distance=0.8;
    end
    
    
    if isempty(criteria.cell_radius_edge)
        criteria.cell_radius_edge=3;
    end
    
    if isempty(criteria.cell_size)
        criteria.cell_size=3; 
    end
    
    
    if isempty(criteria.tumor_cen)
        %criteria.tumor_cen=[117.1875 93.7500 97.5000];
        criteria.tumor_cen=[100.0000 100.000 20.00];

    end
    
    if isempty(criteria.surf_cen)
        %criteria.surf_cen=[110.1174 55.4880 106.7469];
        criteria.surf_cen=[100.0000 100.0000 20.000];
    end
    
    if isempty(criteria.R1)
        criteria.R1=30;
    end
    
    if isempty(criteria.R2)
        criteria.R2=5;
    end

    if isempty(criteria.S1)
        criteria.S1=3;
    end
    
    if isempty(criteria.S2)
        criteria.S2=3;
    end
    
    if isempty(criteria.ref_ratio)
       criteria.ref_ratio=1; 
    end

    criteria
    tic; 
    !rm ~/outmesh.mesh
    [elmtemp,nodtemp]=MakeTetraMesh(eb,pb,criteria);

else 
    tic
    !rm ~/outmesh.mesh
    [elmtemp,nodtemp]=MakeTetraMesh(eb,pb);

end

%% Generating Boundary Nodes and Elements
%Switched the nodes to match the hexahedral meshing!-LS 01/19/16

 nodtemp2(:,2)=nodtemp(:,1)*1e-3; 
 nodtemp2(:,1)=nodtemp(:,2)*1e-3; 
% nodtemp2(:,1)=nodtemp(:,1)*1e-3;
% nodtemp2(:,2)=nodtemp(:,2)*1e-3;

if min(nodtemp(:,3))<0
    nodtemp2(:,3)=nodtemp(:,3)*1e-3+min(nodtemp(:,3))*1e-3+(2/3)*slicethickness*1e-3;
else
nodtemp2(:,3)=nodtemp(:,3)*1e-3+(2/3)*slicethickness*1e-3;
end

disp('Hello');

%[belmtemp bnod]=boundfaces(nodtemp2,elmtemp);
[belmtemp]=getBdyFromMesh(elmtemp(:,1:4)); 
bnod2=unique(belmtemp); 

%%[kk v]=boundary(nodtemp2(:,1),nodtemp2(:,2),nodtemp2(:,3));
%%k=vertcat(kk(:,1),kk(:,2),kk(:,3)); 
%%bnod2=unique(k);
clear bnod;
bnod(:,1:3)=nodtemp2(bnod2,1:3);

%save meshingoutput elmtemp nodtemp2 belmtemp bnod

subplot(1,2,2);
fig.Position=[200 700 1100 1000];
plot3(nodtemp2(:,1),nodtemp2(:,2),nodtemp2(:,3),'.b')
hold on;
plot3(bnod(:,1),bnod(:,2),bnod(:,3),'.r');
title('All Nodes');
xlabel('X'); ylabel('Y'); zlabel('Z')
grid on; AXIS = ([0 max(nodtemp2(:,1)) 0 max(nodtemp2(:,2)) 0 max(nodtemp2(:,3))]);
axis equal; axis(AXIS);
legend('Interior Nodes','Boundary Nodes','Location','Best')


%% Writing everything to files

if ~isempty(elmtemp);
    disp('There are elements!');
    elm=zeros(size(elmtemp,1),6);
    elm(:,1)=1:1:size(elmtemp,1);
    elm(:,2)=elmtemp(:,1); elm(:,3)=elmtemp(:,2); elm(:,4)=elmtemp(:,3); 
    elm(:,5)=elmtemp(:,4); elm(:,6)=elmtemp(:,5);

else 
    elm=elmtemp;
    warning('There are no elements! Remesh using different criteria.');
end

nod=zeros(size(nodtemp2,1),4);
nod(:,1)=1:1:size(nodtemp2,1);
nod(:,2)=nodtemp2(:,1); 
nod(:,3)=nodtemp2(:,2); 
nod(:,4)=nodtemp2(:,3);

% Writing Nodes to File
nodoutf=[outstm '.nod'];
fidnodout=fopen(nodoutf,'w');
fprintf(fidnodout, '%7i %15.8e %15.8e %15.8e 1d0\n',nod');
fclose(fidnodout);

% Writing Elements to File
elmoutf=[outstm,'.elm'];
fidelmout=fopen(elmoutf,'w');
fprintf(fidelmout, '%7i %7i %7i %7i %7i %7i \n',elm');
fclose(fidelmout);

if ~isempty(belmtemp);
    disp('There are boundary elements!');
    belm_ones=(ones(size(belmtemp,1),1))';
    belm_zeros=(zeros(size(belmtemp,1),1))';
    belm_nn=(1:1:size(belmtemp,1))';

    bel(:,1)=belm_nn;
    bel(:,2)=belmtemp(:,1);
    bel(:,3)=belmtemp(:,2);
    bel(:,4)=belmtemp(:,3);
    bel(:,5)=belm_ones;
    bel(:,6)=belm_zeros;

    % Writing Boundary Elements to File
    belf=[outstm '.bel'];
    fidbelf=fopen(belf,'w');
    fprintf(fidbelf,'%7i %7i %7i %7i %7i %7i \n',bel');
    fclose(fidbelf);
else 
    warning('There does not seem to be any boundary elements')
end


if ~isempty(bnod);
    disp('There are boundary nodes!');
    for ii=1:size(bnod,1)
        bnod(ii,4)=find(bnod(ii,1)==nodtemp2(:,1) & bnod(ii,2)==nodtemp2(:,2) & bnod(ii,3)==nodtemp2(:,3));
    end

    bnod_nn=1:1:size(bnod,1);
    bnodfin(:,1)=bnod_nn';
    bnodfin(:,2)=bnod(:,4);

    % Writing Boundary Nodes to File
    bnodf=[outstm,'.bnod'];
    fidbnodf=fopen(bnodf,'w');
    fprintf(fidbnodf,'%7i %7i \n',bnodfin');
    fclose(fidbnodf)
else 
    warning('There does not seem to be any boundary nodes')
    
end




%% Interpolate displacements from grid to the FE mesh
%load MRE_3DMotionData A P MagIm% variables A, P, MagIm
%load HeaderData

A(isnan(A))=0;
P(isnan(P))=0;

fid1=fopen('freq.tmp','wt');

%w=create new file for writing
%t=open file in text mode

fprintf(fid1,'%.4f',freqHz);

%gets 4 places after the decimal point from fid1

fclose(fid1);

%load HeaderData DirIndex %variable DirIndex
%load Mask.mat mask
 
 gen3ddisp(nod, A,P,MREMotionInfo,outstm,mask); %Internal Function
% gen3ddisp(nod, A,P,MagIm,DirIndex,outstm,mask); %Internal Function
% 
 eval(sprintf(['! awk ''{print $1,0,0}'' < ',outstm,'.nod > pressure.s1c']));
% 
% %generating file containing initial guesses
 eval(sprintf(['! awk ''{print $1,1,1,1}'' < ',outstm,'.nod > ',outstm,'.% montagestack(.mtr']));


    % INTERPOLATE DISPLACEMENT DATA ONTO FINITE ELEMENT MESH, GENERATE
    % RECONSTRUCTION ALGORITHIM input FILES, AND LAUNCH MRE
    % RECONSTRUCTION (PARALLEL VERSION)
    
    
    fprintf('Enter Mesh Resolution: \n');
    res = fgetl(fid);
    disp(['Mesh Resolution --> ',res]);
    
    fprintf('\nEnter your email address: \n');
    %read email address
    add = fgetl(fid);
    disp(['Email Address --> ',add]);

    freq = load('freq.tmp');

    fprintf('\nEnter spatial filter weight [0-1]: \n');
    %read spatial filter weight
    spfilt = fgetl(fid);
    disp(['Spatial Filter Weight --> ',spfilt]);

    fprintf('\nEnter maximum number of iterations: \n');
    %read max iterations
    itmax = fgetl(fid);
    disp(['Maximum Iterations --> ',itmax]);

    fprintf('\nEnter regularization parameter: \n');
    %read regularization parameter
    reg = fgetl(fid);
    disp(['Regularization Parameter --> ',reg]);

    fprintf('\nEnter number of processors: \n');
    %read number of processors
    nprocs = fgetl(fid);
    disp(['Number of Processors --> ',nprocs]);


    fprintf('\nPERFORMING DATA INTERPOLATION\n');
    %interpolate data onto finite element mesh using MATLAB scripts

    nod = load([outstm,'.nod']);
    %interpolate measured displacements onto finite element mesh
    gen3ddisp(nod,A,P,MREMotionInfo,outstm,mask);

    %gen3ddisp(nod,A,P,MagIm,DirIndex,outstm,mask);
    
    %generate pressure file
    eval(sprintf(['! awk ''{print $1,0,0}'' < ',outstm,'.nod > pressure.s1c']));

    %generate complex input files  

    cmplx_input(outstm);

    
    %generate file containing initial guess
    eval(sprintf(['! awk ''{print $1,1,1,1}'' < ',outstm,'.nod > ',outstm,'.inv.mtr']));
    
    % POROELASTIC RECONSTRUCTION
    fprintf('\nGENERATING 3D_MRPE_FILES.DAT\n');
   
%     
%     ! mkdir -p MRPE
      ! mkdir -p tet
      
      
      
%     ! mkdir -p MRPE/INV
%     ! mv *.ci pressure.s1c 3D_MRPE* MRPE-* MRPE_*  MRPE/ 
      destination=strcat('tet/',outstm);
      movefile('*.ci',destination);
      movefile('pressure.s1c',destination);
%      ! mv *.ci pressure.s1c tet/

  
    
%% Generating Region Stack File (Adapted from MRIhexmex_v7p3)

%dim=size(MagIm);
xin=0:vox(1):vox(1)*(mridim(1)-1);
yin=0:vox(2):vox(2)*(mridim(2)-1);
zin=0:vox(3):vox(3)*(mridim(3)-1);


   
disp('Mesh Created!');

disp(sprintf('Number of Boundary Nodes: %d \r',size(bnod,1))); 
disp(sprintf('Number of Nodes: %d \r',size(nod,1)));
disp(sprintf('Number of Internal Nodes: %d \r',size(nod,1)-size(bnod,1)))


% Meshing Properties
%[meshprop]=MeshingProp(usedef,default,freqHz,strat);

% output file stem
%[outstm]=Output_File_Sem(usedef,strat,meshprop);

% Displacement Scaling - Many MRE regularization techniques are sensitive
% to the size of the displacements. Either the regularization weights can
% be altered for each case, or the displacements can be scaled so that they
% are always almost the same size.
dispscale = true; % Switch to turn on displacement scaling
dispscalar = 1e-3; % Average displacement amplitude is scaled to this size.
t0=tic;
%% Scale Data if appropriate
meanA = mean(A(repmat(mask,[1 1 1 3])==1));
disp(['Mean Displacement Amplitude all directions ' sprintf('%10.3e',meanA)])
if(dispscale)
    disp(['Scaling Displacements to an average size of ' sprintf('%10.3e',dispscalar) 'm'])
    A=A./meanA.*dispscalar;
end


%% Output mesh

%[namefile]=OutputMesh_v8(outstm,mask,node,deplacement,mridim,coord,regs);

%SaveInputs(outstm,strat,interpo,node,coord,dispscale,dispscalar,default,meshprop,deplacement,msk)

%% Output Runfiles and submit files

% Make sure mu and rho estimates are there
% if(strat.mesh~=2)&&(strat.zone~=2)
%     meshprop.muest=3300;         % Default shear modulus
%     default.rhoest=1000;        % Default density
% end

% Create Directories if they dont exist
direc=pwd;

inpath=[direc '/tet/' outstm '/'];
outpath=[inpath,'inv/'];

% Make hex directory
if (exist('tet','dir')~=7)
    disp('Creating tet directory')
    mkdir('tet')
end

% Make subdirectiories
if (exist(inpath,'dir')~=7)
    disp('Creating tet directory')
    mkdir(inpath)
end
if (exist(outpath,'dir')~=7)
    disp('Creating tet directory')
    mkdir(outpath)
end


% Initialize list of files to move
mvfiles(1).name=[outstm '*']; % Mesh files


%% MREv8 runfile
mudef=3300;
% Estimated properties.    
if(~usedef) 
    muest=input(['Estimate of Shear Modulus (Default ' num2str(mudef) ') >>']);
end
if ~exist('muest','var')||(isempty(muest))
    muest=mudef;
end    
meshprop.muest=muest;

%[zoneprop]=ZoneSize(nod,strat,meshprop);

zoneprop.znedgelength=[1.95633e-02 1.95633e-02 1.95633e-02];
zoneprop.znovlp=[0.15 0.15 0.15];
namefile.regoutf_tet=[outstm '.regstack'];


default.DR=0.18;
default.rhoest=1000;
default.vcomp=0.4;
default.vincomp=0.499;

[regs,noreg]=PriorSegment_v8(MREMotionInfo);  
regoutf=[outstm '.regstack'];
xin=0:vox(1):vox(1)*(mridim(1)-1);
yin=0:vox(2):vox(2)*(mridim(2)-1);
zin=0:vox(3):vox(3)*(mridim(3)-1);

fid=fopen(regoutf,'w');
fprintf(fid,'Array Dimensions \n');
fprintf(fid,'%7i %7i %7i \n',mridim);
fprintf(fid,'x coordinates \n');
fclose(fid);
dlmwrite(regoutf,xin,'-append','delimiter',' ');
fid=fopen(regoutf,'a'); fprintf(fid,'y coordniates \n'); fclose(fid);
dlmwrite(regoutf,yin,'-append','delimiter',' ');
fid=fopen(regoutf,'a'); fprintf(fid,'z coordniates \n'); fclose(fid);
dlmwrite(regoutf,zin,'-append','delimiter',' ');
for ii=1:mridim(3)
    dlmwrite(regoutf,regs(:,:,ii),'-append','delimiter',' ');
end

%% MREv8 runfile

MREMotionInfo(1).tetmesh=1;
MREMotionInfo(1).hexmesh=0;
% Iso incompressible run file - viscoelastic - No soft Prior
%[mvfiles]=RunfileSPoff(outstm,meshprop,mvfiles,namefile,freqHz,outpath,default,vox,zoneprop);
if(~noreg) % Priors is activated. Do not output soft prior runfiles without a supplied segmentation.
    %[mvfiles]=RunfileSPon(outstm,default,meshprop,mvfiles,namefile,freqHz,outpath,vox,zoneprop);
RunfileSPon_v8(MREMotionInfo,outstm,default,meshprop,mvfiles,namefile,freqHz,outpath,vox,zoneprop);
    
end

% Iso incompressible run file - viscoelastic - No soft Prior
% [mvfiles]=RunfileSPoff(outstm,meshprop,mvfiles,namefile,freqHz,outpath,default,vox,zoneprop);
RunfileSPoff_v8(MREMotionInfo,outstm,meshprop,mvfiles,namefile,freqHz,outpath,default,vox,zoneprop);

% Move all files to where they should be.
for ii=1:length(mvfiles)
    [success,message,messageid]=movefile(mvfiles(ii).name,inpath);
    if success==0
        disp(' ');
        disp(['File Transfer Unsuccessful, ii=',num2str(ii)]);
        disp(mvfiles(ii).name);
        disp('MESSAGE: ');
        disp(message);
        disp('MESSAGE ID: ');
        disp(messageid);
    end
end


%tout=toc(touti);
ttotal=toc(t0);

%% Display time for each part:
%disp(['Time for Displacment processing: ' sprintf('%6.2f',tdisp) ' seconds'])
%disp(['Time for FE meshing Process: ' sprintf('%6.2f',temps.tfemesh) ' seconds'])
%disp(['  Time to build incidence list: ' sprintf('%6.2f',tin) ' seconds'])
%disp(['  Time to renumber nodes: ' sprintf('%6.2f',tnodrenum) ' seconds'])
%disp(['  Time to find boundary nodes: ' sprintf('%6.2f',tbnod) ' seconds'])
%disp(['Time to output files: ' sprintf('%6.2f',tout) ' seconds'])
%disp(' ')
disp(['Total processing time: ' sprintf('%6.2f',ttotal) ' seconds'])
disp(['All done with tetrahedral meshing.'])
end




%% Internal Function gen3ddisp   
function gen3ddisp(nod,A,P,MREMotionInfo,outstm,mask);
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   Adapted from MREmesh2.m
%   GEN3DDISP.M
%     This routine is used to compute and interpolate the measured
%     displacement field onto the finite element mesh
%
%   inputS
%     mask        [nx,ny,nz] array of segmented binary images outlining
%                   the region of interest to be meshed
%     stackInfo   a structure containing the pixel-spacing and distance
%                   between slices for the image stack
%     outstm        series outstm
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% pixel spacing
dx=MREMotionInfo(1).pixelspacing_x*1e-3;
dy=MREMotionInfo(1).pixelspacing_y*1e-3;
dz=(MREMotionInfo(1).slicethickness+MREMotionInfo(1).slicegap)*1e-3;

% displacements - stack without top&bottom slices
%Dr=1e-6.*A.*cos(P);
%Di=1e-6.*A.*sin(P);
%A=imrotate(A,-90);
%P=imrotate(P,-90);
A=permute(A,[2 1 3 4]);
P=permute(P,[2 1 3 4]); 
Dr=1e-3.*A.*cos(P);
Di=1e-3.*A.*sin(P);
nx=MREMotionInfo(1).nX;
ny=MREMotionInfo(1).nY; 
nslice=MREMotionInfo(1).nS; 
nd=size(Dr,4); 
clear A P;

% 'pad' displ
pDr=zeros(ny,nx,nslice+2,nd);
pDi=zeros(ny,nx,nslice+2,nd);

pDr(:,:,1,:)=Dr(:,:,1,:);
pDr(:,:,2:nslice+1,:)=Dr;
pDr(:,:,nslice+2,:)=Dr(:,:,end,:);

pDi(:,:,1,:)=Di(:,:,1,:);
pDi(:,:,2:nslice+1,:)=Di;
pDi(:,:,nslice+2,:)=Di(:,:,end,:);

% % Dilate motions around edges a few times to reduce edge interpolation warnings
disp('Dilating to avoid edge warnings')
ndil=3;
mskpad=zeros(ny,nx,nslice+2);
%mskpad(:,:,2:nslice+1,:)=mask;
mskpad(:,:,2:nslice+1,:)=permute(mask,[2 1 3]); 
for ii=1:ndil
    for jj=1:3
        maskdisp=((pDr(:,:,:,jj)~=0)|mskpad);
        [pDr(:,:,:,jj)]=dilatearray3D(pDr(:,:,:,jj),maskdisp);
        [pDi(:,:,:,jj)]=dilatearray3D(pDi(:,:,:,jj),maskdisp);         
    end
end


%[xx,yy,zz]=(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice)*dz);
[xx,yy,zz]=meshgrid(dx:dx:(nx)*dx,dy:dy:(ny)*dy,0:dz:(nslice+1)*dz);

%% Interpolate displacements from grid to the FE mesh

% real displacements
%displ(:,1)=nod(:,1);
displ(:,1)=interp3(xx,yy,zz,pDr(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,2)=interp3(xx,yy,zz,pDr(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
displ(:,3)=interp3(xx,yy,zz,pDr(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

% imaginary displacements
%cdispl(:,1)=nod(:,1);
cdispl(:,1)=interp3(xx,yy,zz,pDi(:,:,:,1),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,2)=interp3(xx,yy,zz,pDi(:,:,:,2),nod(:,2),nod(:,3),nod(:,4),'spline');
cdispl(:,3)=interp3(xx,yy,zz,pDi(:,:,:,3),nod(:,2),nod(:,3),nod(:,4),'spline');

%%Transform Nodal Coordinates according with DirIndex directions
%nod_transform(nod,displ,cdispl,DirIndex,outstm);
nod_transform(nod,displ,cdispl,MREMotionInfo,outstm);

end


%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   CMPLX_input
%     This routine is used to rewrite the complex displacement and
%     pore-pressure data files such that they can be read into FORTRAN
%     as a complex variables
%
%   inputS
%     outstm        series outstm
%
%   Phillip R. Perriñez, Ph.D.
%   Thayer School of Engineering
%   Dartmouth College
%   July 2009
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

function cmplx_input(outstm)

%load .nod file
nod = load([outstm,'.nod']);
nn = size(nod,1);

%load complex displacement file:
UC = load([outstm,'.dsp']);
PC = load('pressure.s1c');

disp_output = [outstm,'.v3c.ci'];
press_output = 'pressure.s1c.ci';

Uxr = UC(:,2);  Uxi = UC(:,3);
Uyr = UC(:,4);  Uyi = UC(:,5);
Uzr = UC(:,6);  Uzi = UC(:,7);

Pr = PC(:,2);   Pi = PC(:,3);

fid=fopen(disp_output,'wt');
d = [Uxr Uxi Uyr Uyi Uzr Uzi];
fprintf(fid,'%6d  (%24.16E,%24.16E)  (%24.16E,%24.16E)  (%24.16E,%24.16E)\n',[1:nn; d']);
fclose(fid);

fid=fopen(press_output,'wt');
p = [Pr Pi];
fprintf(fid,'%6d  (%24.16E,%24.16E)\n',[1:nn; p']);
fclose(fid);

end


%% Nodal Transform
%function nod_transform(nod,displ,cdispl,DirIndex,outstm)
function nod_transform(nod,displ,cdispl,MREMotionInfo,outstm)
%Adapted from MREmesh2    
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
%   NOD_TRANSFORM.M
%     This routine is used to transform the nodal co-ordinate and
%     displacement data using DirIndex
%
%   inputS
%     nod         [nn,4] array containing nodal coordiantes
%     displ       [nn,4] array containing real-valued displacements
%     cdispl      [nn,3] array containing complex-valued displacements
%     DirIndex    4 X 6 matrix containing information regarding how to
%                   transform displacement data and nodal co-ordinates
%     outstm        series outstm
%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo    

% Write Transformed Displacement Vectors to File
fid=fopen([outstm,'.v3r'],'wt');
% fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E\n',new_disp');
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E\n',displ');
fclose(fid);

if MREMotionInfo(1).ExtInt==1 && exist('HeaderData.mat')
    load HeaderData.mat 
    DirIndex=DirIndex(1:3,4:6);
    rowswap=[0 1 0; 1 0 0 ; 0 0 1]; % This makes [1,:,:] = x, [:,1,:]=y, [:,:,1]=z
    MPSto123=rowswap*DirIndex;
    % Perform coordinate transformation
    displ=displ*MPSto123';
    cdispl=cdispl*MPSto123';  % Nnx3*3x3=Nn*3, Transpose IS required
elseif MREMotionInfo(1).ExtInt==1    
    [DirIndex]=DirRead2(MREMotionInfo(1)); 
    rowswap=[0 1 0; 1 0 0 ; 0 0 1]; % This makes [1,:,:] = x, [:,1,:]=y, [:,:,1]=z
    MPSto123=rowswap*DirIndex;
    DirIndex
    % Perform coordinate transformation
    displ=displ*MPSto123';
    cdispl=cdispl*MPSto123';  % Nnx3*3x3=Nn*3, Transpose IS required
end 

%temp = [new_disp(:,1:2) new_cdisp(:,2) new_disp(:,3) new_cdisp(:,3) new_disp(:,4) new_cdisp(:,4)];
temp = [nod(:,1) displ(:,1) cdispl(:,1) displ(:,2) cdispl(:,2) displ(:,3) cdispl(:,3)];

%temp=[displ(:,1) displ(:,3) cdispl(:,2) displ(:,2) cdispl(:,3) displ(:,4) cdispl(:,4)];

% Write Transformed Complex Valued Displacements to File
fid=fopen([outstm,'.dsp'],'wt');
fidv3c=fopen([outstm,'.v3c'],'wt'); %v3c is the same as dsp, just different tags
fprintf(fid,'%6d  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E\n',temp');
fprintf(fidv3c,'%6d  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E  %24.16E\n',temp');
fclose(fid); fclose(fidv3c);

    
end

function [arrayout]=dilatearray3D(arrayin,mask)
% Dilate the edge values in an array to avoid warnings at the boundary during
% interpolation

s=size(arrayin);

arrayout=arrayin;

mask=double(mask);
if(max(mask(:))~=1)
    warning('Dilate3d: Max(mask) ~=1')
end
    

for ii=1:s(1)
    for jj=1:s(2)
        for kk=1:s(3)
            if(mask(ii,jj,kk)~=1)
                msk=mask(max(1,ii-1):min(s(1),ii+1),max(1,jj-1):min(s(2),jj+1),max(1,kk-1):min(s(3),kk+1));
                if(sum(msk(:)>0))
                    ary=arrayin(max(1,ii-1):min(s(1),ii+1),max(1,jj-1):min(s(2),jj+1),max(1,kk-1):min(s(3),kk+1));
                    arrayout(ii,jj,kk)=sum(msk(:).*ary(:))/sum(msk(:));
                end
            end
        end
    end
end

end

function [regs,noreg]=PriorSegment_v8(MREMotionInfo)

if MREMotionInfo(1).SP==1 && isempty(MREMotionInfo(1).SPRegion)==1
    disp('No Soft Prior Regions Defined');
    noreg=1;
    regs=zeros(MREMotionInfo(1).nY, MREMotionInfo(1).nY, MREMotionInfo(1).nS);
else 
    noreg=0;
    spregs=zeros(MREMotionInfo(1).nX, MREMotionInfo(1).nY, MREMotionInfo(1).nS);
    for ii=1:size(MREMotionInfo,2)
        if ~isempty(MREMotionInfo(ii).SPRegion); 
            clearvars regs
            load(MREMotionInfo(ii).SPRegion);
            if exist('regs','var')==0
                error('Define a region for Priors using variable named "regs" with 1s and 0s')
            else 
                spregs(regs(:,:,:)==1)=ii;
            end
        else
            disp(sprintf('You have defined %0d homogenous regions',ii-1));
        end
    end
    clearvars regs
    regs=spregs; 
    montagestack(regs); colorbar; title('SP Regions Defined');
end
end

