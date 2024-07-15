function MRE_recon_3D_2014_08_07_v08(recon_type,nc_rank)
% INPUTS:
%   folder - folder path to put data.
%   data_type - what type of data: '2D', '3D', and '3D_nav'
%   recon_type - type of recon to run

% Anh: 09/12/10
% Use min-var phase error correction
% Support flexible echoes and spiral in/out readouts
% Joe:2011_01_05 make changes in order to automate recon
% Joe:2011_06_07 Convert to full brain with multislab to match sequence.
%                Only handles single echo.
% Joe:2011_10_06 Convert to use GPU, must use FM and Sense
% Joe:2011_11_10 Shift Image to center of the field of view before recon
% Joe:2011_11_10 Add recon type feature
% Joe:2011_12_06 Navigator Recon can use SENSE
% Joe:2011_12_06 Enable Multiple GPUs
% Joe:2012_06_05 Enable reconstruction of additional image sizes
% Joe:2012_07_11 Use full pathnames for loading SENSE and FM maps

%suppress warnings for printing to screen 
%#ok<*PRTCAL>
if ~(isempty(dir('*.dat')))
    fname_d = dir('*.dat');
    fname = fname_d(1).name;
    fname = fname(1:end-4);
else
     sprintf('Did not find a data file to reconstruct \n')
     return
end



fprintf(1, 'Starting Image Reconstruction\n')
tic

switch recon_type
    case 1 % typical recon
        perform_motion_correction = 1;
        perform_gridding_recon = 0;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;
        CPU_par = 1;
        short_ntp = 0;
        skip_ntp = 0;
    case 2 % fast recon
        perform_motion_correction = 1;
        perform_gridding_recon = 0;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;   
        CPU_par = 1;
        short_ntp = 1;
        skip_ntp = 0;
    case 3 % gridding recon
        perform_motion_correction = 0;
        perform_gridding_recon = 1;
        perform_iterative_recon = 0;
        perform_field_correction = 0;
        mask_recon = 0;
        sense_recon = 0;   
        CPU_par = 0;
        short_ntp = 0;
        skip_ntp = 0;
    case 4
        perform_motion_correction = 1;
        perform_gridding_recon = 0;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;   
        CPU_par = 1;
        short_ntp = 0;
        skip_ntp = 1;
    case 5
        perform_motion_correction = 1;
        perform_gridding_recon = 1;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;
        CPU_par = 1;
        short_ntp = 0;
        skip_ntp = 0;
    case 6
        perform_motion_correction =0;      
        perform_gridding_recon = 0;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;
        CPU_par = 1;
        short_ntp = 0;
        skip_ntp = 0;
    case 7 % everything
        perform_motion_correction = 1;
        perform_gridding_recon = 0;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;
        CPU_par = 1;
        short_ntp = 0;
        skip_ntp = 0;
    otherwise
        fprintf(1, 'No valid recon type specified\n')
        perform_motion_correction =     1;
        perform_gridding_recon = 1;
        perform_iterative_recon = 1;
        perform_field_correction = 1;
        mask_recon = 1;
        sense_recon = 1;
        CPU_par = 0;
        short_ntp = 0;
        skip_ntp = 0;
end

    nav_on = 1;
    bSingle_shot = 0;


if CPU_par
    fprintf(1, '   Parallel CPUs will be used\n')
end
if perform_gridding_recon
    fprintf(1, '   Gridding reconstruction will be run\n')
end
if perform_iterative_recon
    fprintf(1, '   Iterative reconstruction will be run\n')
end
if perform_motion_correction
    if nav_on
        fprintf(1, '    -Motion correction willl be used\n')
    else
        fprintf(1, '    -WARNING: NO NAVIGATOR DATA\n')
        fprintf(1, '       Motion Corrrection has been turned off\n')
        perform_motion_correction = 0;
    end
end
if perform_field_correction
    fprintf(1, '    -field correction will be used\n')
end
if sense_recon
    fprintf(1, '    -a SENSE reconstruction will be performed\n')
end


if CPU_par == 1
    if matlabpool('size') == 0
        matlabpool open
    end
% elseif par > 1
%     if ~(matlabpool('size') == par)
%         if ~(matlabpool('size') == 0)
%             matlabpool close;
%         end
%         matlabpool(par)
%     end
else
    if ~(matlabpool('size') == 0)
        matlabpool close;
    end
end

folder = strcat(pwd,'/');
fname = strcat(folder,fname,'.dat');

if exist('ascconv.mat','file')
    load ascconv.mat
else
    parsasc_full(fname);
    load ascconv.mat
end


ntp = ascconv.lRepetitions +1;
if short_ntp
    ntp_recon =3;
else
    ntp_recon=ntp;
end

if skip_ntp
    ntp_start =20;
else
    ntp_start=1;
end

nc = size(ascconv.asCoilSelectMeas.asList,2);
D = ascconv.sSliceArray.asSlice(1).dReadoutFOV/10;
N = ascconv.sKSpace.lBaseResolution;
nl = ascconv.sKSpace.lPhaseEncodingLines;
nl_design = ascconv.sWiPMemBlock.adFree(1);
% nl_design = 4; %THIS IS FROM SEQUENCE BUG -- NEED TO PUT BACK IN WIPMEMBLOCK
nslab = ascconv.sGroupArray.asGroup.nSize;
niter = 20;

dummies = 1;


if ~exist('nc_rank','var')
    nc_rank = nc;
end

mkdir recon_info
[read_shift, phase_shift] = image_shift(ascconv);
save recon_info/corners.mat read_shift phase_shift

% read_shift = 0;


    nsl = ascconv.sKSpace.lPartitions;
    df =  -ascconv.sGroupArray.asGroup.dDistFact;  
    % df = -1/4;
% df =1/3;
%     df = 0;

N_recon = N;
nsl_recon = nsl;

% exc_order = [1 3 5 7 2 4 6 8];
exc_order = [0 2 4 6 8 1 3 5 7 9];
% exc_order = [0 2 4 6 8 10 1 3 5 7 9 11];
% exc_order = [0 2 8 4 6 10 1 3 9 5 7 11];
%exc_order = [0 2 4 6 8 10 12 14 1 3 5 7 9 11 13 15];

if ~exist(strcat(folder,'recon_info'),'dir')
    mkdir recon_info
end
if ~exist(strcat(folder,'raw_data'),'dir')
    mkdir raw_data
end

%% K-space trajectory for high res data.
if exist(strcat(folder,'recon_info/kspace_info.mat'),'file')
    load recon_info/kspace_info.mat
    fprintf(1,'K-space information loaded\n')
else
    fprintf(1,'Setting up K-space information\n')
    gamp = 2.1; %   was always 2.4
    gslew = 130; %200;  % was 120
    gts = 10e-06;
    tsamp = 5e-06;
    fprintf(1,' -Calculating k-space trajectory\n')
    [Gx,Gy,kxi,kyi,sx,sy] = genspi(D,N,nl_design,gamp,gslew,gts);
    kxt=interp1([0:gts:gts*length(kxi)-gts],kxi,[0:tsamp:gts*length(kxi)-tsamp])';
    kyt=interp1([0:gts:gts*length(kyi)-gts],kyi,[0:tsamp:gts*length(kyi)-tsamp])';

    nk = length(kxt)-2;
    kx = zeros(nk,nl);
    ky = zeros(nk,nl);
    kxo = kxt(1:nk);
    kyo = kyt(1:nk);
    ng = length(Gx);
    gx = zeros(ng,nl);
    gy = zeros(ng,nl);
    %rotate matrix for proper orientation
    sprintf('Performing %d rotations',nl)
    phi = 2*pi/nl;
%     phi = 2*pi/nl_design;
    for ii = 0:(nl-1)
         ang_rot = phi*(ii-(nl-1)*floor(ii/nl));
         kx(:,ii+1) = kxo*cos(ang_rot) + kyo*sin(ang_rot);
         ky(:,ii+1) = kyo*cos(ang_rot) - kxo*sin(ang_rot);
         gx(:,ii+1) = Gx*cos(ang_rot) + Gy*sin(ang_rot);
         gy(:,ii+1) = Gy*cos(ang_rot) - Gx*sin(ang_rot);
    end
    kx = kx(:);
    ky = ky(:);
    gx = gx(:);
    gy = gy(:);
    kx = reshape(kx,length(kx(:))/nl,nl); 
    ky = reshape(ky,length(ky(:))/nl,nl); 
    kz =  -(-floor(nsl/2):1:(ceil(nsl/2) - 1));
    ww =  weight_vor(col(kx),col(ky),nl,1);    % decimated DCF
    thresh=find(ww>1, 1 );
    for inl=1:nl
        ww((inl-1)*length(ww)/nl+thresh:inl*length(ww)/nl)=1;
    end
    nro = length(kx(:))/nl;
    Nro_adc = 512; 
    N_adc = ceil(nro/Nro_adc);
    % Two times oversampling in readout direction, real and imaginary
    lgadc = ceil(length(gx)/nl*(gts/tsamp));
    nroa = 8*ceil(lgadc/8);
    
    % READOUT PHASE
    gamma = 2*pi*4.257e3;
    gambar = gamma/(2*pi);
    rosc = (gamp*D*gambar*tsamp);
    rophs = repmat(col(0:(length(kx)-1))*rosc,[1 2]);

    if nav_on
        % K-space trajectory for low res data.
        fprintf(1,' -Calculating navigator k-space trajectory\n')
%         LR_SIZE_Z = 10;
%         LR_SIZE_Z = 6;
        LR_SIZE_Z = 8;
        LR_SIZE_XY = 15; %8
        LR_nl = 1;      % number of interleave for low res.
        [kxLRd, kyLRd, gxLRd, gyLRd]=genkspace_mod(D,LR_SIZE_XY,0,LR_nl,gamp,gslew,tsamp,0,0,gts,0,1,4,floor(1),1); %#ok<NASGU>
        kxLRd = reshape(kxLRd,length(kxLRd(:))/LR_nl,LR_nl); 
        kyLRd = reshape(kyLRd,length(kyLRd(:))/LR_nl,LR_nl);
        gxLRd = reshape(gxLRd,length(gxLRd(:))/LR_nl,LR_nl); 
        kxLR = kxLRd(:, 1);
        kyLR = kyLRd(:, 1);
        gxLR = gxLRd(:, 1);
        wwLR = weight_vor(kxLR(:),kyLR(:),1,1);    % decimated DCF
        wwLR(wwLR>1)=1;
        wwLR(isnan(wwLR))=1;
%         kzLR = -(-5:1:4);       % check with pulse sequence!
        kzLR = -(-4:1:3);       % check with pulse sequence!
        echoZ = floor(nsl/2) + 1;
        start_idx = - floor(LR_SIZE_Z/2) + echoZ;
        stop_idx = ceil(LR_SIZE_Z/2) - 1 + echoZ;
        % LowRes
        nroLR = length(kxLR(:))/1;
        lgadcLR = ceil(length(gxLR)/1*(gts/tsamp));
        nroaLR = 8*ceil(lgadcLR/8);
        
        % READOUT PHASE
        roscLR = (gamp*D*gambar*tsamp);
        rophsLR = col(0:(length(kxLR)-1))*roscLR;
        
        save recon_info/kspace_info kx ky kz gx kxLR kyLR kzLR gxLR ww wwLR start_idx stop_idx nro gts LR_SIZE_XY LR_SIZE_Z tsamp N_adc Nro_adc nl nl_design lgadc nroa nroLR lgadcLR nroaLR rophs rophsLR
    else
        save recon_info/kspace_info kx ky kz gx ww nro gts tsamp N_adc Nro_adc nl nl_design lgadc nroa rophs
    end  
    fprintf(1,' -K-space information saved\n')  
end

%% Read in data
%do not read data if data has already beeen read in
if ~exist(strcat(folder,'raw_data/raw_data_1_1.mat'),'file')
    fprintf(1,'Reading in Raw Data\n')

    fid = fopen(fname,'r');

    bytes_to_skip_hdr = fread(fid,1,'uint32');
    fseek(fid,bytes_to_skip_hdr,'bof');

    curvol_tmp = zeros(nro, nl, nsl, nslab, nc);
    if nav_on
        nav_vol_tmp = zeros(nroLR, LR_SIZE_Z, nl,nc,nsl,nslab);
    end
    for kk = 1: ntp
        fprintf(1,' -Reading data for image %d\n',kk)
        for jj = 1:nsl
            for ll = 1:nl
                for ii = 1:nslab 
                    for dd = 1:N_adc
                        for cc = 1:nc
                            % Read the data from here.
                            sMDH = ice_read_mdh_va21(fid);
                            if ~(sMDH.ushSamplesInScan == Nro_adc)
                                sprintf('HELP: WRONG NUMBER OF SAMPLES')
                                keyboard
                            end
                            [raw,rcount] = fread(fid,2*Nro_adc,'float32'); %#ok<NASGU>
                            if dd < N_adc
                                curvol_tmp(((dd-1)*Nro_adc+1):(dd*Nro_adc),ll,jj,ii, cc) = raw(1:2:2*Nro_adc)+1i*raw(2:2:2*Nro_adc);
                            else
                                n_extra = mod(nro,Nro_adc);
                                if n_extra == 0
                                    n_extra = Nro_adc;
                                end
                                curvol_tmp(((dd-1)*Nro_adc+1):(((dd-1)*Nro_adc+1)+n_extra-1),ll,jj,ii, cc) = raw(1:2:2*n_extra)+1i*raw(2:2:2*n_extra);
                            end
                        end
                    end
                    if nav_on
                        for idxn = 1:LR_SIZE_Z     % Read the navigator data here.
                            for cc = 1:nc
                                sMDH = ice_read_mdh_va21(fid);     
                                %check navigator read in data!!!!!!!!!!
                                if ~(sMDH.ushSamplesInScan == nroaLR)
                                    sprintf('HELP: WRONG NUMBER OF SAMPLES IN NAVIGATOR!');
                                    keyboard
                                end
                                [raw, rcount] = fread(fid, 2*nroaLR, 'float32'); %#ok<NASGU>
                                nav_vol_tmp(:, idxn, ll,cc,jj,ii) = raw(1:2:2*nroLR) + 1i*raw(2:2:2*nroLR);
                            end
                        end
                    end
                end
            end
        end
        fprintf(1,sprintf(' -Saving data %d\n',kk))
        for ii = 1:nslab
            curvol = permute(curvol_tmp(:,:,:,ii,:),[1 2 3 5 4]);
            eval(sprintf('save raw_data/raw_data_%d_%d curvol',ii,kk)) 
            if nav_on
                nav_vol = nav_vol_tmp(:,:,:,:,:,ii);
                eval(sprintf('save raw_data/raw_nav_%d_%d nav_vol',ii,kk))
            end
        end
    end
    fclose(fid);
    clear curvol sMDH bytes_to_skip_hdr curvol_tmp
    if nav_on
        clear nav_vol nav_vol_tmp
    end
    fprintf(1,' -Data read complete\n')
else
    %fprintf(1,'Data loaded\n')
end



%% Gridding recon
if perform_gridding_recon
    %FIXME: are the gridding results shifted by one slice?
    fprintf(1,'Performing gridding reconstruction\n')
    img = zeros(N_recon,N_recon,nsl,nc,nslab);
    img_grid = zeros(N_recon,N_recon,nsl*nslab);
    for kk=1:ntp
         for ii = 1:nslab
           eval(sprintf('load raw_data/raw_data_%d_%d',ii,kk));  
%          for ii =3
             fprintf(1,' -Gridding recon for slab %d of %d image %d of %d starting\n',ii,nslab,kk,ntp)
             if nsl > 1
                 for cc = 1:nc
                     t_vol = fftshift(ifft(fftshift(curvol(:,:,:, cc), 3), [], 3), 3);
                     for jj = 1:nsl
                         img(:,:, jj, cc, ii) = k2image(col(kx), col(ky), col(t_vol(:, :, jj).*exp(-1i*(read_shift.*kx+phase_shift.*ky)*2*pi)), col(ww), N_recon, 2.5); 
                     end
                 end
             else
                 parfor cc = 1:nc
                    img(:,:, 1, cc, ii) = k2image(col(kx), col(ky), col(curvol(:,:,1, cc).*exp(-1i*(read_shift.*kx+phase_shift.*ky)*2*pi)), col(ww), N_recon, 2.5); 
                 end
             end
             img_grid(:,:, (exc_order(ii) )*nsl+1: (exc_order(ii)+1)*nsl)=flipdim(sqrt(sum(abs(img(:,:,:,:, ii)).^2, 4)),3);
         end
         eval(sprintf('save img_grid_%d img_grid',kk))
    end
    clear img img_grid t_vol curvol
    fprintf(1,'     gridding recon complete\n')
end

%% Set up SENSE and FM maps
if sense_recon||perform_field_correction
    fprintf(1, 'Setting up SENSE map and field map\n')
    if ((~exist(strcat(folder,'SE_FM/FM.mat'),'file'))||(~exist(strcat(folder,'SE_FM/sen.mat'),'file')))
        fprintf(1, ' -creating SENSE map and field map\n')
        senmap_run3D('SE_FM',read_shift,phase_shift)
    end
    nsl_overlap = round(nslab*nsl*(1-df)+(nsl-nsl*(1-df)));
    nsl_overlap_recon = round(nslab*nsl_recon*(1-df)+(nsl_recon-nsl_recon*(1-df)));
    if perform_field_correction
        if ~exist(strcat(folder,'recon_info/FM_new_1.mat'),'file')
            fprintf(1, ' -Resampling field map\n')
            if ~exist('FM','var')
                eval(sprintf('load %sSE_FM/FM',folder))
            end
            if 0
                FM_tmp = zeros(size(FM,1),size(FM,1),2*size(FM,3));
                FM_tmp(:,:,17:48) = FM;
                FM = FM_tmp;
            end
            FM = resample_map_resolution(FM,N,nsl_overlap); %#ok<NODEF>
            FM_tmp = FM;
            for ii = 1:nslab
                start_idx = exc_order(ii)*nsl*(1-df)+1;
                stop_idx = start_idx + nsl -1;
                FM = FM_tmp(:,:,start_idx:stop_idx);
                eval(sprintf('save recon_info/FM_new_%d FM',ii))
            end
            fprintf(1, ' -New field maps saved\n')
        end
        if (N ~= N_recon)||(nsl~=nsl_recon)
            if ~exist(strcat(folder,'recon_info/FM_new_1_',num2str(N_recon),'_',num2str(N_recon),'_',num2str(nsl_recon),'.mat'),'file')
                eval(sprintf('load %sSE_FM/FM',folder))
                FM = resample_map_resolution(FM,N_recon,nsl_overlap_recon);
                FM_tmp = FM;
                for ii = 1:nslab
                    start_idx = exc_order(ii)*nsl_recon*(1-df)+1;
                    stop_idx = start_idx + nsl_recon -1;
                    FM = FM_tmp(:,:,start_idx:stop_idx);
                    eval(sprintf('save recon_info/FM_new_%d_%d_%d_%d FM',ii,N_recon,N_recon,nsl_recon))
                end
                fprintf(1, ' -New field maps saved\n')
            end
        end
    end
    clear FM FM_tmp
    if mask_recon
        if ~exist(strcat(folder,'recon_info/mask_new_1.mat'),'file')
            fprintf(1, ' -Resampling mask\n')
            if ~exist('mask','var')
                eval(sprintf('load %sSE_FM/mask',folder))
            end
            if 0
                mask_tmp = zeros(size(mask,1),size(mask,1),2*size(mask,3));
                mask_tmp(:,:,17:48) = mask;
                mask = mask_tmp;
            end
            mask = resample_map_resolution(mask,N,nsl_overlap);
            mask_tmp = mask;
            for ii = 1:nslab
                start_idx = exc_order(ii)*nsl*(1-df)+1;
                stop_idx = start_idx + nsl -1;
                mask = mask_tmp(:,:,start_idx:stop_idx);
                eval(sprintf('save recon_info/mask_new_%d mask',ii))
            end
            fprintf(1, ' -New masks saved\n')
        end
        if (N ~= N_recon)||(nsl~=nsl_recon)
            if ~exist(strcat(folder,'recon_info/mask_new_1_',num2str(N_recon),'_',num2str(N_recon),'_',num2str(nsl_recon),'.mat'),'file')
                eval(sprintf('load %sSE_FM/mask',folder))
                mask = resample_map_resolution(mask,N_recon,nsl_overlap_recon);
                mask_tmp = mask;
                for ii = 1:nslab
                    start_idx = exc_order(ii)*nsl_recon*(1-df)+1;
                    stop_idx = start_idx + nsl_recon -1;
                    mask = mask_tmp(:,:,start_idx:stop_idx);
                    eval(sprintf('save recon_info/mask_new_%d_%d_%d_%d mask',ii,N_recon,N_recon,nsl_recon))
                end
                fprintf(1, ' -New masks saved\n')
            end
        end
    end
    clear mask mask_tmp
    if sense_recon
        if ~exist(strcat(folder,'recon_info/sen_new_1.mat'),'file')
            fprintf(1, ' -Resampling SENSE map\n')
            if ~exist('sen','var')
                eval(sprintf('load %sSE_FM/sen',folder))
            end
            if 0
                sen_tmp = ones(size(sen,1),size(sen,2),2*size(sen,3),nc);
                sen_tmp(:,:,17:48,:) = sen;
                sen = sen_tmp;
            end
            sen_tmp = zeros(N,N,nsl_overlap,nc);
            for cc = 1:nc
                sen_tmp(:,:,:,cc) = resample_map_resolution(sen(:,:,:,cc),N,nsl_overlap); %#ok<NODEF>
            end
            for ii = 1:nslab
                start_idx = exc_order(ii)*nsl*(1-df)+1;
                stop_idx = start_idx + nsl -1;
                sen = sen_tmp(:,:,start_idx:stop_idx,:);
%                 sensz = size(senx);
%                 [~,s,v] = svd(reshape(senx,[],nc),0);
%                 sen = reshape(reshape(senx,[],nc)*v(:,1:nc_rank),[sensz(1) sensz(2) sensz(3) nc_rank]);
                eval(sprintf('save recon_info/sen_new_%d sen',ii))
            end
            fprintf(1, ' -New SENSE maps saved\n')
        end
        if (N ~= N_recon)||(nsl~=nsl_recon)
             if ~exist(strcat(folder,'recon_info/sen_new_1_',num2str(N_recon),'_',num2str(N_recon),'_',num2str(nsl_recon),'.mat'),'file')
                 eval(sprintf('load %sSE_FM/sen',folder))
                 sen_tmp = zeros(N_recon,N_recon,nsl_overlap_recon,nc);
                for cc = 1:nc
                    sen_tmp(:,:,:,cc) = resample_map_resolution(sen(:,:,:,cc),N_recon,nsl_overlap_recon);
                end
                for ii = 1:nslab
                    start_idx = exc_order(ii)*nsl_recon*(1-df)+1;
                    stop_idx = start_idx + nsl_recon -1;
                    sen = sen_tmp(:,:,start_idx:stop_idx,:);
                    eval(sprintf('save recon_info/sen_new_%d_%d_%d_%d sen',ii,N_recon,N_recon,nsl_recon))
                end
                fprintf(1, ' -New SENSE maps saved\n')
             end
        end
    end
    clear sen sen_tmp start_idx stop_idx nsl_overlap nsl_overlap_recon
end

%% Perform motion correction
if perform_motion_correction
    if ~exist(strcat(folder,'recon_info/kspace_Corrections.mat'),'file')
        moco3D_meanshot
        !cp recon_info/kspace_Corrections_meanshot.mat recon_info/kspace_Corrections.mat
    else
        fprintf(1,'Motion correction has previously been run\n')
    end
end

%% Iterative recon for DW images
if perform_iterative_recon

    iterativestart=toc;
    fprintf(1,'Beginning Iterative Reconstruction\n')
    % Cn=0;
    beta = 1e-8;
    delta = 1e-7;
    % Cn = beta*C3D_sparse(ones(N,N,nsl));

    if perform_motion_correction
        fprintf(1,' -Loading motion corrections\n')
        load recon_info/kspace_Corrections
    else
        corrections = zeros(4,nl, nsl, nslab, ntp);  
    end   
    if perform_field_correction
        L = ceil(nro/500);
    else
        L = 0;
    end
    tt = ((0:(nro-1))*tsamp)';
    tt_new = repmat(tt,[1 nl nsl]);
    mkx = zeros(nro,nl,nsl);
    mky = zeros(nro,nl,nsl);
    mkz = zeros(nro,nl,nsl);
    imginitA = zeros(N_recon, N_recon, nsl_recon); 
    % mask = col(ones(size(imginit)));
    cimg = zeros(N_recon, N_recon, nsl_recon*nslab);
    for jj = 1:nsl
        for ll = 1:nl
            mkz(:,ll,jj) = kz(jj);
            mkx(:,ll,jj) = kx(:,ll);
            mky(:,ll,jj) = ky(:,ll);
        end
    end

    
    sen_tmp_all = zeros(N_recon*N_recon*nsl_recon,nc,nslab);
    FM_tmp_all = zeros(N_recon*N_recon*nsl_recon,nslab);
    mask_tmp_all = zeros(N_recon*N_recon*nsl_recon,nslab);
    for ii = 1:nslab
        sen_tmp_all(:,:,ii) = sense_map( ii, N, nsl,nc, folder,N_recon,nsl_recon);
        FM_tmp_all(:,ii) = field_correction(perform_field_correction, ii, N, nsl, folder,N_recon,nsl_recon);
        mask_tmp_all(:,ii) = data_mask(mask_recon, ii, N, nsl, folder,N_recon,nsl_recon);
    end
    
        c_tmp = zeros(N_recon,N_recon,nsl_recon,nslab);
        for kk = (ntp_start+dummies):ntp_recon
                parfor ii=1:nslab
                    sen_tmp = sen_tmp_all(:,:,ii);
                    FM_tmp = FM_tmp_all(:,ii);
                    mask_tmp = mask_tmp_all(:,ii);
                    
                    imginit = col(imginitA);
                    
                    R = Robject(mask_tmp,'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');%,'huber','delta',delta);
                    
                    fprintf(1,' -Iterative recon for slab %d of %d image %d of %d starting\n',ii,nslab,kk,ntp)
                    mkx_tmp = zeros(nro,nl,nsl);
                    mky_tmp = zeros(nro,nl,nsl);
                    mkz_tmp = zeros(nro,nl,nsl);
                    curvol = get_raw_data(ii,kk);
                    for cc = 1:nc
                        curvol(:,:,:,cc) = curvol(:,:,:,cc).*exp(-1i*(read_shift.*mkx+phase_shift.*mky)*2*pi);  %shift image in plane and by half a pixel in the z direction
                        % curvol(:,:,:,cc) = curvol(:,:,:,cc).*exp(-1i*(read_shift.*mkx+phase_shift.*mky+1/(2*nsl).*mkz)*2*pi);  %shift image in plane and by half a pixel in the z direction

                    end

                    for ll = 1:nl
                        for jj =1:nsl
                            mkz_tmp(:,ll,jj) = mkz(:,ll,jj) +corrections(3,ll,jj,ii,kk) ;
                            mkx_tmp(:,ll,jj) = mkx(:,ll,jj)+corrections(1,ll,jj,ii,kk) ;
                            mky_tmp(:,ll,jj) = mky(:,ll,jj)+corrections(2,ll,jj,ii,kk);
                            curvol(:,ll,jj,:) = curvol(:,ll,jj,:).*exp(1i*corrections(4,ll,jj,ii,kk));  
                        end
                    end
                    
                    A = fast_mr_v2(col(mkx_tmp), col(mky_tmp), col(mkz_tmp), D, N_recon, N_recon, nsl_recon,2*N_recon, 2*N_recon, 2*nsl_recon, 5, col(tt_new), FM_tmp(logical(mask_tmp)), 0, L, 1,0,0,logical(mask_tmp));     
                    % A = fast_mr_3D(col(mkx_tmp), col(mky_tmp), col(mkz_tmp), D, N_recon, N_recon, nsl_recon,nl,nsl_recon, 2*N_recon, 2*N_recon, 2*nsl_recon, 5, tt, FM_tmp_all(:,ii), 0, L, 1);     
                    S = sense_svd(A ,sen_tmp(logical(mask_tmp),:),nc_rank); % S = sense(A ,sen_map);
                    t_img = pwls_pcg1(imginit(logical(mask_tmp)), S, 1,prepData(S,col(curvol)), R, 'niter', niter);
                    % t_img = qpwls_pcg_Anh(imginit(logical(mask_tmp)), S, 1,prepData(S,col(curvol)), 0, Cn, 1, niter); 
                    c_tmp(:,:,:,ii) = reshape(embed(t_img,logical(mask_tmp)), N_recon, N_recon, nsl_recon);
                end

            % end 
            for ii=1:nslab
                cimg(:,:, (exc_order(ii) )*nsl_recon+1: (exc_order(ii)+1)*nsl_recon)= c_tmp(:,:,:,ii);
            end
            eval(sprintf('save img_%d  cimg', kk));
            
%             start_idx = round(nsl*(-df)/2);
%             stop_idx = nsl - start_idx;
%             slices_keep = stop_idx-start_idx;
%             for ii = 1:nslab
%                 img(:,:,((ii-1)*slices_keep+1):(slices_keep*ii)) = cimg(:,:,((ii-1)*nsl+start_idx+1):((ii-1)*nsl+stop_idx));
%             end
%             eval(sprintf('save img_cut_%i img',kk))

            MRE_3D_cut_single(cimg,kk)
            
        end
    
    fprintf(1,' -iterative recon complete\n')
    % MRE_3D_cut(ntp)
end

%% Cleanup code
time_end = toc;
fprintf(1,'Image reconstruction complete\n')
fprintf(1,'Total time = %d\n',time_end)
if perform_iterative_recon
    fprintf(1,'iterative time = %d\n',time_end-iterativestart)
end

clear all;
end

%% functions
function we = field_correction(bOn, slab_num, N, nsl, folder, N_recon, nsl_recon)
    if bOn
        if (N_recon~=N)||(nsl~=nsl_recon)
            eval(sprintf('load %srecon_info/FM_new_%d_%d_%d_%d',folder,slab_num,N_recon,N_recon,nsl_recon));
        else
            eval(sprintf('load %srecon_info/FM_new_%d',folder,slab_num));
        end
        we = col(FM);
    else
        we = col(zeros(N_recon, N_recon, nsl_recon));
    end
    clear FM
    
end

function we = field_correction_nav(bOn, slab_num, N, nsl, folder)
    if bOn
        eval(sprintf('load %srecon_info/FM_new_nav_%d',folder,slab_num));
        we = col(FM);
    else
        we = col(zeros(N, N, nsl));
    end
    clear FM
end

function msk = data_mask(msk_on, slab_num, N, nsl, folder, N_recon, nsl_recon)
    if msk_on
        if (N_recon~=N)||(nsl~=nsl_recon)
            eval(sprintf('load %srecon_info/mask_new_%d_%d_%d_%d',folder,slab_num,N_recon,N_recon,nsl_recon));
        else
            eval(sprintf('load %srecon_info/mask_new_%d',folder,slab_num));
        end
        msk = col(mask);
    else
        msk = col(zeros(N_recon, N_recon, nsl_recon));
    end
    clear mask
    
end

function sen_map = sense_map( slab_num, N, nsl,nc, folder, N_recon, nsl_recon)
    if (N_recon~=N)||(nsl~=nsl_recon)
        eval(sprintf('load %srecon_info/sen_new_%d_%d_%d_%d',folder,slab_num,N_recon,N_recon,nsl_recon));
    else
        eval(sprintf('load %srecon_info/sen_new_%d',folder,slab_num));
    end
    sen_map = reshape(sen,N_recon*N_recon*nsl_recon,nc);   
    sen_map(find(abs(sen_map)>2)) = 2;
    clear sen
end

function sen_map = sense_map_nav( slab_num, N, nsl,nc, folder)
    eval(sprintf('load %srecon_info/sen_new_nav_%d',folder,slab_num));
    sen_map = reshape(sen,N*N*nsl,nc);   
    sen_map(find(abs(sen_map)>2)) =2;
    clear sen
end

function raw = get_raw_data(slab_num, image_num)
    eval(sprintf('load raw_data/raw_data_%d_%d',slab_num,image_num));
    raw = curvol;
    clear curvol
end

function [free_GPU_index status img] = get_free_GPU_id(GPU_id_list,GPUs_in_use,index, Nx,Ny,Nz)
    data_folder = '/data/jholtrop/GPU/tmp_GPU';
    out_file = 'out.file';
    get_data = 0;
    status = 0; %returns 1 when a completed job was removed and will be returned
    img = 0;
    free_GPU_index = 0; %negative means there are no free GPUs, otherwise a valid GPU id is returned
    if index > 0 %check if a specific GPU id has finished
        if exist(sprintf('%s%i/%s',data_folder,GPU_id_list(index),out_file),'file')
            get_data =1;
            GPU_data_to_get = GPU_id_list(index);
            free_GPU_index = index; 
        elseif GPUs_in_use(index) == 0;
            free_GPU_index = index; 
        end    
    else   %Want to find a free GPU ID
        for ii=1:length(GPU_id_list)
            if exist(sprintf('%s%i/%s',data_folder,GPU_id_list(ii),out_file),'file')
                get_data =1;
                GPU_data_to_get = GPU_id_list(ii);
                free_GPU_index = ii;
                break
            elseif GPUs_in_use(ii) == 0;
                free_GPU_index = ii;
                break
            end    
        end
        if free_GPU_index == 0
            pause(5); %wait 5 seconds
        end
    end
    if get_data
        img = read_3D_GPU_multi(Nx,Ny,Nz,GPU_data_to_get);
        status = 1;
    end
end