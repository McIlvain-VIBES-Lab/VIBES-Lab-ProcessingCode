function [propstack_RE,propstack_IM] = forward_mtr_files_to_MR_interp(dspf,dim,res)
%forward_mtr_files_to_MR_interp Interpolates the forward material files back to an MRE stack.
%   Inputs: dspf: Name of displacement file, this function will interoplaote the
%   associated .forward.mtr files. The property mesh is what the properties
%   are actually stored on for the forward and inverse problem. 
%   There was a slight bug in the forward problem of MREv9p21a and below 
%   which leaves the forward .meshind file empty. It is fixed as of 3/14/19. 
%   If the .meshind file is empty, you will need to edit it and fill it in. 
%   dim and res <[x y z]> the matrix size and voxel size (m) of the stack
%   to interpolate to. If they are not supplied, an MRE file (NLI_Input.mat
%   or MRE_3DMotionData.mat) file can be selected from the user input. 
%   Outputs: mat file with same filestem as the dsp file.
%            propstack_RE and propstack_IM: interpolated properties real
%            and imag parts, propstack(:,:,:
%   I also fixed a bug in v9p33, which did not output the nodesense arrays
%   so cant see where the boundaries of the displacement mesh are.  
if(nargin<1)
    [dspf,dsppath]=uigetfile('*.dsp','Select the hexahedral displacement file');    
else
    dsppath=pwd;
end

if(nargin<2)
    [MRf,MRpath]=uigetfile('*','Select MRE file to interpolate back to');
    load(fullfile(MRpath,MRf));
    if(exist('A','var')) % Old MRE_3DMotionData format
        dim=[size(A,1) size(A,2) size(A,3)];
        load(fullfile(MRpath,'HeaderData.mat'));
        res=DirIndex(4,1:3)/1000;
    elseif(exist('Ur','var'))
        dim=[size(Ur,1) size(Ur,2) size(Ur,3)];
        res=voxsize_mm/1000;        
    end
end


% Find all the property files
fstm=fullfile(dsppath,dspf);
fstm=fstm(1:end-12); % This excludes the '.' 

meshind=load([fstm '.meshind']);
if(isempty(meshind))
    error([fstm '.meshind is empty. MREv9p21a and above fixes this as of 3/13/09. Otherwise fill in the meshind file'])
end

nprops=size(meshind,2); 

propstack_RE=zeros([dim nprops]);
propstack_IM=zeros([dim nprops]);

for ii=1:nprops
    % Interpolate real properties
    mshf=[fstm '.mtrmesh.' sprintf('%2.2i',meshind(1,ii)) '.nod'];
    mtrf=[fstm '.RE..forward.prop.' sprintf('%2.2i',ii) '.mtr'];
    disp(['Prop ' int2str(ii) ' real: '])
    disp(['File name : ' mtrf])
    disp(['mesh name : ' mshf])
    [propstack_RE(:,:,:,ii)]=hexfile_to_MR_interp(mtrf,mshf,dim,res);
    
    % Interpolate Imag properties
    mshf=[fstm '.mtrmesh.' sprintf('%2.2i',meshind(2,ii)) '.nod'];
    mtrf=[fstm '.IM..forward.prop.' sprintf('%2.2i',ii) '.mtr'];
    disp(['Prop ' int2str(ii) ' imag: '])
    disp(['File name : ' mtrf])
    disp(['mesh name : ' mshf])
    [propstack_IM(:,:,:,ii)]=hexfile_to_MR_interp(mtrf,mshf,dim,res);
end
save([fullfile(dsppath,dspf(1:end-3)) 'props.mat'],'propstack_RE','propstack_IM','dim','res') 


end


