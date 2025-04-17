%*********************************************************************
%   tet_output_general.m
%   Modified from shear_output by Matt McGarry
%   Original file written by
%   AJ Pattison    
%
%*********************************************************************
%   function stacks=tet_output(filename,nodf,showfig,S,stack)
%   Inputs: filename: Columns of this file will be interploated.
%           showfig: (optional, default 'y') - if 'y',  If 'clip', the stack is clipped so 
%                    that zero portions are not plotted. Otherwise no 
%                    values are plotted. all 5 arguments must be supplied 
%                    if showfig is specified.
%           S: (optional, default = [] (auto)) - Montage Dimensions.
%           stack: (optional, default 0) size of each image slice. Set 
%                  stack = 0 to use the same resolution as the Motion data, 
%                  other values will use interpolation.
%
%   Output: stacks: Structure, length = number of data columns (node number
%                   column ignored). stacks(1).vals = interpolated first
%                   data column, stacks(2).vals = interpolated 2nd data
%                   column, etc.....
%
%   Examples: stacks = tet_output('filename.v3c') 
%               will plot the values in the tetrahedral motion data file 
%               'filename.v3c', and save them as a structure, stacks. 
%               Default montage size = auto calculated. Default 
%               value of stack=0, therefore the resolution is the same as 
%               the original data. 
%             stacks = tet_output('INV/reconresult.99','y',[3 3]) 
%               will plot the values in the tetrahedral reconstructed file 
%               'INV/reconresult.99', and save them as a structure, stacks.
%               Resolution will be same as data, the montage will be 3x3.
%             stacks = tet_output('dataset.v3c','clip',[3 3]) 
%               will plot the values in the tetrahedral motion data file 
%               'dataset.v3c', and save them as a structure, stacks. 
%               Montage will be 3x3, Images will be trimmed to only show 
%               non-zero data.
%             stacks = tet_output('dataset.v3c','n',[3 3],500) 
%               will not plot the values, only save them as a structure, 
%               stacks. stack=500, therefore reolution will be interpolated
%               to 500.
%
%   the variable stacks is also saved in a mat-file 'tetresults.mat'
%
%   This function reformats the input data and calls the C mex to 
%   create a shear image stack based on the shear values of the nodes.
%
%   NOTE: the mex file requires a lot of memory especially when the 
%   mesh is highly refined, check your compiler for memory capability.
%
%   To compile the .mex file: actualize the path to the mex file
%   location and enter 'mex -o shear_stack_mex.c' in the MATLAB command
%   window.
%
%
%********************************************************************

function stacks=tet_output_general(filename,nodf,elmf,cols,stack_dim,voxel_size,showfig,S)


if(nargin==0)
    [filename,pathname]=uigetfile('*','Select tetrahedral value file to interpolate');
    filename=fullfile(pathname,filename);
end
if(nargin<=1) % Get name of nod file
    [nodf,pathname]=uigetfile('*.nod','Select appropriate node file');
    nodf=fullfile(pathname,nodf);
end
if(nargin<=2) % Get name of elm file
    [elmf,pathname]=uigetfile('*.elm','Select appropriate element file');
    elmf=fullfile(pathname,elmf);
end
% Check size of input files
vals=load(filename);
nod=load(nodf);
nod = nod(:,2:4)*1000; % converting nodal positions to mm
elm=load(elmf);
elm = elm(:,2:5); % element incidence list

if(size(vals,1)~=size(nod,1))
    disp(['Tet value file: ' int2str(size(vals,1)) ' Nodes, tet node file ' int2str(size(nod,1)) ' Nodes']);
    error('Tet value file must have same number of entries as tet nod file')
end

if(nargin<=3)
    cols=2:size(vals,2);
end

if(nargin<=4) % Find stack_dim and voxel_size from MRE data
    [mref,pathname]=uigetfile('*.mat','Select MRE data file to interpolate back to');
    load(fullfile(pathname,mref));
    if(exist('A','var')) % MRE_3DMotionData format
        stack_dim=[size(A,1) siz(A,2) size(A,3)];
        [headerf,pathname]=uigetfile('*.mat','Select HeaderData file with resolution info');
        load(fullfile(pathname,headerf));
        voxel_size=DirIndex(4,1:3);
    elseif(exist('Ur','var')) % NLI standard format
        stack_dim=[size(Ur,1) size(Ur,2) size(Ur,3)];
        voxel_size=voxsize_mm;
    end
end
    
if(nargin<7)
    showfig='y';
elseif(nargin<8)
    S=[];
end

stack_size = zeros(stack_dim);

%% call mex file
nv=0;
for ii=cols % 1:ncols
    tetstack = shear_stack_mex(stack_size, voxel_size, int32(elm), nod, vals(:,ii));
    nv=nv+1;
    stacks(nv).vals=tetstack;
    %% display
    if (strcmp(showfig,'y'))          
        montagestack(stacks(nv).vals);drawnow;pause(0.2);
        title(['\bfvalue in column ',int2str(ii)])
        colorbar
        drawnow;pause(0.1);
    end
end
save tetresults stacks
end

