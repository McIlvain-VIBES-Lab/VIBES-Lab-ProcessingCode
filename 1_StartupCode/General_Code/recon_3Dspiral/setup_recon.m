function recoInfo = setup_recon(fname,traj_type,rec_type, gamp, gslew)%, IPAT_fact, shft_FOVx, shft_FOVy, C_regularize)
% function recoInfo = setup_recon(fname,traj_type,rec_type)
%
%  setup reconstruction: Can use [] as any input to get default behavior
%
%  Inputs
%        fname - filename of data, such as meas.dat. 
%		Can also be [] and will find the .dat file or .out 
%		file in the current directory
%	 traj_type: 	1: constant density spiral
%			2: Variable density, original (Kim, 2xos design)
%			3: Variable density, modified
%	 rec_type:	1: gridding, individual coil images (k2image)
%			2: gridding, sum-of-squares         (k2image)
%			3: gridding with field map (FM), individual coils	(k2image_we)
%			4: gridding with field map (FM), sum-of-squares        	(k2image_we)
%			5: iterative with field map (FM), individual coils     	(fast_mr)
%			6: iterative with field map (FM), sum-of-squares       	(fast_mr)
%			7: iterative with SENSE, 		(fast_mr)
%			8: iterative with SENSE and FM		(fast_mr)		
%			9: (future) iterative with FM gradients
%			10: (future) SPIRIT (non-Cartesian Grappa)
%                       11: (future) 3D
%                       12: (future) 3D with 3D navigator
%                       13: Motion induced phase cor (DTI) with 2D
%                       14: Motion induced phase cor (DTI) with 3D
%        (optional)
%	  gamp     :   maximum gradient amplitude for spiral design (default 2.2 G/cm = 22 mT/m)
%         glsew    :   maximum gradient slew rate for spiral design (default 140 mT/m/ms)
%         IPAT_fact  : reduction factor for reduced FOV applications, can be non-integer if traj_type is 3.
%         shft_FOVx  : number of pixels to shift the FOV by in the x-direction
%         shft_FOVy  : number of pixels to shift the FOV by in the y-direction
%         C_regularize : C penalty matrix for iterative reconstruction with quadratic penalty
%
%  Outputs
%	  recoInfo - a structure containing information needed for reconstruction
%
%
%
%
%       WARNINGS: 
%       1. reduced FOV acquisitions with iterative reconstruction are poorly conditioned. 
%		Need to specify region, different recon code required.
%
%
% 	FUTURE TODO:
%       0. include shots to use for sense reconstruction (recon every shot vs every image)
%	1. Get IPAT_fact recon working for all three trajectory types
%	2. handle indexing recon with slice and time point vector, so as to not need to reconstruct entire time series, slices
%	3. make varargin handling, instead of list
%       4. Include eddy current delays
%       5. Fessler field map negative - yep, it is. Has positive sign convention.

% GLOBAL VARIABLES NOT PASSED IN YET
gts = 10e-6;   %Gradient raster time is 10 us
tsamp = 5e-6;  % Sampling rate for ADC is 5 us



% get filename of data file if not specified
flag_findname = 0;
if exist('fname','var')
   if isempty(fname)
      flag_findname = 1;
   end
else
   flag_findname = 1;
end


if flag_findname
   if (exist('meas.dat','file') | exist('meas.out','file'))
        fname = 'meas';
   elseif ~(isempty(dir('*.dat')))
        fname_d = dir('*.dat');
        fname = fname_d(1).name;
        fname = fname(1:end-4);
    elseif ~(isempty(dir('*.out')))
         fname_d = dir('*.out');
         fname = fname_d(1).name;
         fname = fname(1:end-4);
    else
         sprintf('Did not find a data file to reconstruct \n')
         return
    end
end

% Get fname_root and fname for data
    fname_list = dir(sprintf('%s*',fname));
    ii = 1;
    while (fname_list(ii).isdir)
        ii = ii+1;
    end    
    if (~isempty(findstr(fname_list(ii).name,'.out')) | ~isempty(findstr(fname_list(ii).name,'.asc')))
        sw_version = 'VA';
        fname_root = fname;
        fname_data = sprintf('%s.out',fname);
        fname_hdr = sprintf('%s.asc',fname);
    elseif (~isempty(findstr(fname_list(ii).name,'.dat')))
        sw_version = 'VB';
        fname_root = fname;
        fname_data = sprintf('%s.dat',fname);
        fname_hdr = sprintf('%s.dat',fname);
    else
        sprintf('Did not find data file and file name root \n')
        return
    end

%GET actual size of adc
fid = fopen(fname_data);
bytes_to_skip_hdr = fread(fid,1,'uint32');
fseek(fid,bytes_to_skip_hdr,'bof');
sMDH = ice_read_mdh_va21(fid);
nroa = sMDH.ushSamplesInScan; 
fclose(fid);


% set trajectory type if not specified
if ~exist('traj_type','var')
   traj_type = 1; % Default: constant density spiral
elseif isempty(traj_type)
   traj_type = 1; % Default: constant density spiral
end
% set reconstruction type if not specified
if ~exist('rec_type','var') 
   rec_type = 1; % Default: Gridding individual coil images
elseif isempty(rec_type)
   rec_type = 1; % Default: Gridding individual coil images
end
% gradient amplitude
if ~exist('gamp','var') 
   gamp = 2.2; % Default: 2.2 G/cm = 22 mT/m
elseif isempty(gamp)
   gamp = 2.2; % Default: 2.2 G/cm = 22 mT/m
end
% gradient slew rate
if ~exist('gslew','var') 
   gslew = 140; % Default: 140 mT/m/ms
elseif isempty(gslew)
   gslew = 140; % Default: 140 mT/m/ms
end
% % reduced FOV factor
% if ~exist('IPAT_fact','var') 
%    IPAT_fact = 1; % Default: no reduced FOV
% elseif isempty(IPAT_fact)
%    IPAT_fact = 1; % Default: no reduced FOV
% end
% % number of pixels to shift FOV in x
% if ~exist('shft_FOVx','var') 
%    shft_FOVx = 0; % Default: no shifted FOV
% elseif isempty(shft_FOVx)
%    shft_FOVx = 0; % Default: no shifted FOV
% end
% % number of pixels to shift FOV in y
% if ~exist('shft_FOVy','var') 
%    shft_FOVy = 0; % Default: no shifted FOV
% elseif isempty(shft_FOVy)
%    shft_FOVy = 0; % Default: no shifted FOV
% end
% % Regularizer for iterative reconstruction
% if ~exist('C_regularize','var') 
%    C_regularize = 0; % Default: no regularization
% elseif isempty(C_regularize)
%    C_regularize = 0; % Default: no regularization
% end

nsl = 1; % To initialize in case slices not defined in header
      


% Load sequence parameters, run parsasc if needed (ie no seqinfo.m)
if exist('ascconv.mat','file')
    load ascconv.mat
else
    parsasc_full(fname_data);
    load ascconv.mat
end

fov = ascconv.sSliceArray.asSlice(1).dReadoutFOV;
nx = ascconv.sKSpace.lBaseResolution;
ny = ascconv.sKSpace.lPhaseEncodingLines;
nsl = ascconv.sSliceArray.lSize;

ucMode = ascconv.sSliceArray.ucMode;

if ucMode == 1
    slorder = 1:1:nsl;
elseif ucMode == 2
    slorder = nsl:-1:1;
elseif ucMode == 4
    sl1 = flipdim(nsl:-2:1,2);
    sl2 = flipdim((nsl-1):-2:1,2);
    slorder = [sl1 sl2];
else
    slorder = 1:1:nsl;
end
    

slthick = ascconv.sSliceArray.asSlice(1).dThickness;

if isfield(ascconv,'lRepetitions')
    ntp = ascconv.lRepetitions +1;
else
    ntp = 1;
end
num_coils = size(ascconv.asCoilSelectMeas.asList,2);

if (ascconv.asCoilSelectMeas.iUsedRFactor) == 3
    triple_mode = 1;
else
    triple_mode = 0;
end

[shft_FOVx, shft_FOVy] = image_shift(ascconv);

TR = ascconv.alTR;
TE = ascconv.alTE;

FAdeg = ascconv.adFlipAngleDegree;

IPAT_fact = 1;
C_regularize = 0;

% if ~exist('ntp','var')
%   ntp = 1;
% end

% if ~exist('num_coils','var')
%     num_coils=1;
% end
% 
% if ~exist('triple_mode','var')
%     triple_mode=0;
% end


%size of voxels for analyze header
szx = fov/nx;
szy = fov/nx;
szz = slthick;

%NUMBER OF ADCS
% This will be used for long readouts that create multiple ADCs
% num_adc = 1;


%-------------------------------
% CALCULATE k-SPACE TRAJECTORY
%-------------------------------
    D = fov/10;
    N = nx; 
    nl = ny;

     
    switch traj_type
	case 1    % 1: constant density spiral
    	  % Generate k-space design
    	  [kx, ky, gx, gy]=genkspace2(D,N,0,IPAT_fact*nl,gamp,gslew,tsamp,0,0,gts);

          %Accomodate IPAT_fact
          kx = reshape(kx,length(kx(:))/(IPAT_fact*nl),IPAT_fact*nl);
          ky = reshape(ky,length(ky(:))/(IPAT_fact*nl),IPAT_fact*nl);

          kx = kx(:,1:IPAT_fact:end);
          ky = ky(:,1:IPAT_fact:end);
       
    
    	  % Fix array size
    	  kx = col(kx);
    	  ky = col(ky);
     
	  % Two times oversampling in readout direction, real and imaginary
	  nro = length(kx)/nl;
	  lgadc = ceil(length(gx)/(IPAT_fact*nl)*(gts/tsamp));
	  %nroa = 8*ceil(lgadc/8);
           %if (lgadc > 4096)
              num_adc = ceil(lgadc/1024);
           %end

	case 2   %2: Variable density, original (Kim, 2xos design)
            % Generate k-space design
    	  [kx, ky, gx, gy]=genkspace2(D,N,0,IPAT_fact*nl,gamp,gslew,tsamp,0,0,gts,1);

          %Accomodate IPAT_fact
          kx = reshape(kx,length(kx(:))/(IPAT_fact*nl),IPAT_fact*nl);
          ky = reshape(ky,length(ky(:))/(IPAT_fact*nl),IPAT_fact*nl);

          kx = kx(:,1:IPAT_fact:end);
          ky = ky(:,1:IPAT_fact:end);
       
    
    	  % Fix array size
    	  kx = col(kx);
    	  ky = col(ky);
     
	  % Two times oversampling in readout direction, real and imaginary
	  nro = length(kx)/nl;
	  lgadc = ceil(length(gx)/(IPAT_fact*nl)*(gts/tsamp));
	  %nroa = 8*ceil(lgadc/8);
           %if (lgadc > 4096)
              num_adc = ceil(lgadc/1024);
           %end

   	case 3	   %3: Variable density, modified
            sprintf('Traj Type 3 not completed \n')
            return
            % NEED TO ASSIGN: kx,ky, nro, lgadc, (if not 1, then num_adc)

        otherwise
             sprintf('Trajectory Type unknown \n')
             return
         end        
 
	  % Calculate weighting based on reduced set of shots
    	  ww = weight_vor(col(kx),col(ky),nl,1);
    	  ww(find(ww>(3*IPAT_fact)))=0;
    	  ww(find(isnan(ww)))=0;



%---------
% Set up timing vector
%--------------
tt = [];
for ii = 1:nl
    tt = [tt; (0:nro-1).'*tsamp+(TE*1e-6)];
end


%-----------------------------
% NOW FILL IN recoInfo
%-----------------------------

recoInfo.fname_root = fname_root;
recoInfo.fname_data = fname_data;
recoInfo.fname_hdr = fname_hdr;
recoInfo.kx = kx;
recoInfo.ky = ky;
recoInfo.ww = ww;
recoInfo.traj_type = traj_type;
recoInfo.rec_type = rec_type;
recoInfo.shft_FOVx = shft_FOVx;
recoInfo.shft_FOVy = shft_FOVy;
recoInfo.N = N;
recoInfo.FOV = fov;
recoInfo.szx = szx;
recoInfo.szy = szy;
recoInfo.szz = szz;
recoInfo.TR = TR*1e-6;
recoInfo.TE = TE*1e-6;
recoInfo.tt = tt;
recoInfo.FAdeg = FAdeg;
recoInfo.nl = nl;
recoInfo.ntp = ntp;
recoInfo.slthick = slthick;
recoInfo.nsl = nsl;
recoInfo.slorder = slorder;
recoInfo.IPAT_fact = IPAT_fact;
recoInfo.num_coils = num_coils;
recoInfo.triple = triple_mode;
recoInfo.gam = gamp;
recoInfo.gslew = gslew;
recoInfo.gts = gts;
recoInfo.tsamp = tsamp;
recoInfo.nro = nro;
recoInfo.nroa = nroa;
recoInfo.lgadc = lgadc;
recoInfo.num_adc = num_adc;
recoInfo.num_iters = 20;
recoInfo.L = 7;
recoInfo.J = 6;
recoInfo.C_regularize = C_regularize;
recoInfo.shots_to_use = nl; % For sense reconstruction, can use different number of shots for recon
recoInfo.moco = 0;

save recoInfo recoInfo

    
