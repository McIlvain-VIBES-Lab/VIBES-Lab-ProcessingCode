function status = spiral_recon(recoInfo)
%function status = spiralraw(recoInfo)
%   recoInfo from setup_recon.m

recdirlist = {'recon_coils';'recon_sos';'recon_coils_fmcp';'recon_sos_fmcp';'recon_coils_fmiter';'recon_sos_fmiter';'recon_sen';'recon_fm'};
  

sense_reco_list = [7;8];
FM_reco_list = [3;4;5;6;8];

flag_sense = max(ismember(recoInfo.rec_type,sense_reco_list));
flag_FM = max(ismember(recoInfo.rec_type,FM_reco_list));

if flag_sense
   if exist('sen.mat','file')
       load sen
   elseif exist('../sen.mat','file')
       load ../sen
   else
       sprintf('Need sen.mat in top level directory for sense reconstruction \n')
       return
   end

  if ~(size(sen,1) == recoInfo.N)
     %RESAMPLE SEN for matrix size of current reconstruction
      for coilIndex = 1:recoInfo.num_coils
          sen_new(:,:,:,coilIndex) = resample_map_resolution(sen(:,:,:,coilIndex),recoInfo.N);   
      end
      sen=sen_new;
   end

end

if flag_FM
   if exist('FM.mat','file')
      load FM
   elseif exist('../FM.mat','file')
      load ../FM
   else
      sprintf('Need FM.mat in top level directory for field corrected reconstruction \n')
      return
   end
   if ~(size(FM,1) == recoInfo.N)
       %RESAMPLE FM FOR matrix size of current reconstruction
       FM_new = resample_map_resolution(FM,recoInfo.N);
       FM = FM_new;
   end
 
   FM_orig = FM;
end


%keyboard
kx = reshape(recoInfo.kx,[recoInfo.nro recoInfo.nl]);
ky = reshape(recoInfo.ky,[recoInfo.nro recoInfo.nl]);
nro = recoInfo.nro;
nroa = recoInfo.nroa;


skp = 128;


fid = fopen(recoInfo.fname_data);
bytes_to_skip_hdr = fread(fid,1,'uint32');
fseek(fid,bytes_to_skip_hdr,'bof');

bytesPerImage = recoInfo.num_coils*(((skp/2) * 2) + ((2*recoInfo.nroa) * 4)); 

%%%checking file size
%%fsz = (nroa*2*4+128)*nl*nsl*ntp*num_echo;

tmpvol = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl);
dat = zeros(recoInfo.nroa*recoInfo.num_adc,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils,recoInfo.ntp); 
dat = zeros(recoInfo.nro,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils,recoInfo.ntp); 

slorder = recoInfo.slorder;

for imageIndex = 1:(recoInfo.ntp)
    for shotIndex = 1:recoInfo.nl %% load each interleaf's data into the buffer
        for sliceIndex = 1:recoInfo.nsl
            for adcIndex = 1:recoInfo.num_adc
                for coilIndex = 1:recoInfo.num_coils %% read in each coil data
                    
                    %% Read and parse interleaf points
                    sMDH = ice_read_mdh_va21(fid);
                    [raw,rcount] = fread(fid,2*nroa,'float32');
                    
                    if (rcount*recoInfo.num_adc < 2*nro)
                        sprintf('Reached end of file before running out of indices \n')
                        return
                    end
                    
                    %% Put it in the correct buffer position
                    tmp = raw(1:2:2*nroa)+i*raw(2:2:2*nroa);
                    start_ind = (adcIndex-1)*nroa + 1;
                    end_ind = min(adcIndex*nroa,nro);
                    num_ind = end_ind - start_ind +1;
                    tmp2 = col(tmp(1:num_ind)).*exp(-i*2*pi*ky(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVy).*exp(-i*2*pi*kx(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVx);
                    dat(start_ind:end_ind,shotIndex,slorder(sliceIndex),coilIndex,imageIndex) = tmp2;
                    
                end % end coilIndex
            end  % end adcIndex
            
            
            if (recoInfo.num_coils == 12)&&(recoInfo.triple==1)
                dat_tmp = squeeze(dat(:,shotIndex,slorder(sliceIndex),:,imageIndex));
                dat_tmp2 = tim_12ch(dat_tmp);
                for coilIndex = 1:recoInfo.num_coils
                    dat(:,shotIndex,slorder(sliceIndex),coilIndex,imageIndex) = dat_tmp2(:,coilIndex);
                end
            end
            
            
        end  % end sliceIndex
    end  % end shotIndex
end % end imageIndex

corrections = zeros(3,recoInfo.nl, recoInfo.nsl, recoInfo.ntp);
if recoInfo.moco == 1
    if ~exist('corrections.mat','file')
        corrections = mle_moco(recoInfo,dat);
    else
        load corrections.mat
    end
end

if isempty(gcp('nocreate'))
   parpool
end

for reconIndex = 1:length(recoInfo.rec_type)

    recnum = recoInfo.rec_type(reconIndex);
    recondir = recdirlist{recnum};
    
    mkdir(recondir)
    cd(recondir)
    
    save recoInfo.mat recoInfo

% Call correct recon and save output

for imageIndex = 1:(recoInfo.ntp)
    for sliceIndex = 1:recoInfo.nsl
        for shotIndex = 1:recoInfo.nl
            mkx(:,shotIndex,sliceIndex) = kx(:,shotIndex)-corrections(1,shotIndex,sliceIndex,imageIndex);
            mky(:,shotIndex,sliceIndex) = ky(:,shotIndex)-corrections(2,shotIndex,sliceIndex,imageIndex);
            mdat(:,shotIndex,sliceIndex,:,imageIndex) = dat(:,shotIndex,sliceIndex,:,imageIndex).*exp(1i*(-1)*corrections(3,shotIndex,sliceIndex,imageIndex));
        end
    end
    
    switch recnum
        %--------------------------------------------------------------
        case {1,2} %gridding, individual coil images or sos (k2image)
        %--------------------------------------------------------------
            for sliceIndex = 1:recoInfo.nsl
                parfor coilIndex = 1:recoInfo.num_coils
                    imgtmp(:,:,sliceIndex,coilIndex) = k2image(col(mkx(:,:,sliceIndex)),col(mky(:,:,sliceIndex)),col(mdat(:,:,sliceIndex,coilIndex,imageIndex)),col(recoInfo.ww),recoInfo.N,2.5); %uncorrected for inhomogeneity
                end
            end
            
            if (recnum == 1)
                for coilIndex = 1:recoInfo.num_coils
                    img = squeeze(imgtmp(:,:,:,coilIndex));
                    save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
                end
            else
                img = sos_image(imgtmp);
                save(sprintf('sos_%05d',imageIndex),'img')
            end
            
            
            
        %--------------------------------------------------------------
        case {3,4} %gridding with field map (FM), individual coils or sos	(k2image_we)
        %--------------------------------------------------------------
            for sliceIndex = 1:recoInfo.nsl
                parfor coilIndex = 1:recoInfo.num_coils
                    imgtmp(:,:,sliceIndex,coilIndex) = k2image_we(col(mkx(:,:,sliceIndex)),col(mky(:,:,sliceIndex)),col(mdat(:,:,sliceIndex,coilIndex,imageIndex)),col(recoInfo.ww),recoInfo.N,col(FM_orig(:,:,sliceIndex)),recoInfo.L, col(recoInfo.tt),2.5); %corrected for inhomogeneity
                end
            end
            
            if (recnum == 3)
                for coilIndex = 1:recoInfo.num_coils
                    img = squeeze(imgtmp(:,:,:,coilIndex));
                    save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
                end
            else
                img = sos_image(imgtmp);
                save(sprintf('sos_%05d',imageIndex),'img')
            end
            
        %--------------------------------------------------------------
        case {5,6} %iterative with field map (FM), individual coils  or sos  	(fast_mr)
        %--------------------------------------------------------------
            
            
            xinit = col(zeros(recoInfo.N));
            mask = ones(recoInfo.N);
            for sliceIndex = 1:recoInfo.nsl
                A = fast_mr(col(mkx(:,:,sliceIndex)),col(mky(:,:,sliceIndex)), recoInfo.FOV/10, recoInfo.N, 2*recoInfo.N, recoInfo.J, col(recoInfo.tt), col(FM_orig(:,:,sliceIndex)), 0, recoInfo.L, 1);
                
                parfor coilIndex = 1:recoInfo.num_coils
                    recon = qpwls_pcg(xinit, A, 1, col(mdat(:,:,sliceIndex,coilIndex,imageIndex)), 0, recoInfo.C_regularize, 1, recoInfo.num_iters, mask, 0);
                    imgtmp(:,:,sliceIndex,coilIndex) = reshape(recon(:,end), recoInfo.N, recoInfo.N);
                end
            end
            
            if (recnum == 5)
                for coilIndex = 1:recoInfo.num_coils
                    img = squeeze(imgtmp(:,:,:,coilIndex));
                    save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
                end
            else
                img = sos_image(imgtmp);
                save(sprintf('sos_%05d',imageIndex),'img')
            end
            
        %--------------------------------------------------------------
        case {7,8} %iterative with SENSE, and SENSE and FM		(fast_mr)
        %--------------------------------------------------------------
            if (recnum == 7)
                FM = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl);
            else
                FM = FM_orig;
            end
            tt = reshape(recoInfo.tt,[length(recoInfo.tt)/recoInfo.nl recoInfo.nl]);
            
            sen = reshape(sen, [recoInfo.N*recoInfo.N recoInfo.nsl recoInfo.num_coils]);
            
            xinit = col(zeros(recoInfo.N));
            mask = ones(recoInfo.N);
            parfor sliceIndex = 1:recoInfo.nsl
                    A = fast_mr(col(mkx(:,:,sliceIndex)), col(mky(:,:,sliceIndex)), recoInfo.FOV/10, recoInfo.N, 2*recoInfo.N, recoInfo.J, col(tt(:,:)), col(FM(:,:,sliceIndex)), 0, recoInfo.L, 1);
                    As = sense(A,squeeze(sen(:,sliceIndex,:)));
                    recon = qpwls_pcg(xinit, As, 1, col(squeeze(mdat(:,:,sliceIndex,:,imageIndex))), 0, recoInfo.C_regularize, 1, recoInfo.num_iters, mask, 0);
                    imgtmp(:,:,sliceIndex) = reshape(recon(:,end), recoInfo.N, recoInfo.N);
            end
            img = imgtmp;
            save(sprintf('sen_%05d',imageIndex),'img')
            
            
            
            
        %--------------------------------------------------------------
        otherwise
            sprintf('Recon Type %d not implemented yet \n',recoInfo.rec_type)
            return
    end
    
    
    
    
end % end imageIndex

cd ..

end

fclose('all')



